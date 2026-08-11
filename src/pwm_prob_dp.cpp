#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cmath>

// This function computes P(m): the probability that a PWM-width seq window scores
// strictly above the cutoff, given each seq nt (0th-order) or di-nt (1st-order) context.
// PWM scores are discretized onto an integer grid controlled by 'resolution', and
// the summed score distribution is built by convolution (pos by pos), instead of enumerating all 4^motif_width nt paths.
// At each position, shift and weight the existing score probs for A/C/G/T, then combine branches that fall in the same score bin.
// Seqs are processed one at a time, so memory is O(resolution) and independent of nr. seqs.

static constexpr double GRID_EPS = 1e-9; // avoid floating point issues when assigning bin
static const int INTERRUPT_CHECK_INTERVAL = 128; // check for user interrupts every 128 seqs

static inline double floor_grid(double x) {
  return std::floor(x + GRID_EPS); // avoid floating-point error leading to prior bin placement
}

// constrain probabilities to [0, 1] avoiding floating-point rounding errors
static inline double clamp01(double x) {
  return (x < 0.0) ? 0.0 : ((x > 1.0) ? 1.0 : x);
}

// Add the current score probs after shifting them by the nt score and weighting them by the nucleotide probability.
// Combine probability mass for all passing branches in the final score bin.
static inline void add_weighted_score_prob(const double* score_prob,
                                           int n_score_bins, int nt_score_shift, double nt_prob,
                                           double* updated_score_prob, int n_updated_bins) {

  int n_direct_bins = std::min(n_score_bins, n_updated_bins - nt_score_shift);

  if (n_direct_bins < 0) n_direct_bins = 0;

  for (int score_bin = 0; score_bin < n_direct_bins; ++score_bin) {
    updated_score_prob[score_bin + nt_score_shift] += score_prob[score_bin] * nt_prob;
  }

  if (n_direct_bins < n_score_bins) {
    double cutoff_prob = 0.0;

    for (int score_bin = n_direct_bins; score_bin < n_score_bins; ++score_bin) {
      cutoff_prob += score_prob[score_bin];
    }
    updated_score_prob[n_updated_bins - 1] += cutoff_prob * nt_prob;
  }
}

extern "C" SEXP pwm_prob_dp(SEXP pwm_sexp, SEXP cutoff_sexp,
                           SEXP seqs_freq_nt_sexp, SEXP cond_prob_sexp,
                           SEXP resolution_sexp) {
  // R object -> C++ object
  const double* pwm = REAL(pwm_sexp); // pwm: 4 x motif_width, rows A,C,G,T
  const int motif_width = Rf_ncols(pwm_sexp);
  const double cutoff = REAL(cutoff_sexp)[0];
  const double* seqs_freq_nt = REAL(seqs_freq_nt_sexp); // seqs_freq_nt: n_seqs_eval x 4
  const int n_seqs_eval = Rf_nrows(seqs_freq_nt_sexp); // number of seqs being evaluated
  const int resolution = INTEGER(resolution_sexp)[0]; // positive integer
  const bool order1 = (cond_prob_sexp != R_NilValue);
  const double* cond_prob = order1 ? REAL(cond_prob_sexp) : nullptr; // cond_prob: 4 x 4 x n_seqs_eval, or R_NilValue for order 0

  SEXP result = PROTECT(Rf_allocVector(REALSXP, n_seqs_eval));
  double* window_prob = REAL(result); // P(m) per seq

  // use mem managed by R (R_alloc) so it is safely released if the fun call is interrupted
  double* min_score_by_pos = (double*) R_alloc(motif_width, sizeof(double)); // min PWM score per pos
  double* max_score_by_pos = (double*) R_alloc(motif_width, sizeof(double)); // max PWM score per pos
  // find the min and max nt score at each PWM position
  for (int mot_pos = 0; mot_pos < motif_width; ++mot_pos) {
    double min_pos_score = pwm[0 + 4 * mot_pos];
    double max_pos_score = min_pos_score;
    for (int nt = 1; nt < 4; ++nt) {
      const double nt_score = pwm[nt + 4 * mot_pos];
      if (nt_score < min_pos_score) min_pos_score = nt_score;
      if (nt_score > max_pos_score) max_pos_score = nt_score;
    }
    min_score_by_pos[mot_pos] = min_pos_score;
    max_score_by_pos[mot_pos] = max_pos_score;
  }

  // sum the position-specific bounds to get the full-window PWM score bounds
  long double min_pwm_score_sum = 0.0L, max_pwm_score_sum = 0.0L; // long double to match R's sum()
  for (int mot_pos = 0; mot_pos < motif_width; ++mot_pos) {
    min_pwm_score_sum += min_score_by_pos[mot_pos];
    max_pwm_score_sum += max_score_by_pos[mot_pos];
  }
  const double min_pwm_score = static_cast<double>(min_pwm_score_sum);
  const double max_pwm_score = static_cast<double>(max_pwm_score_sum);
  const double score_range = max_pwm_score - min_pwm_score;

  // handle fixed motif; when every window scores identically regardless of nt content.
  if (score_range <= 0.0) {
    // P(m) is 1 if the fixed window score passes the strict cutoff, otherwise 0.
    const double deterministic_window_prob = (min_pwm_score > cutoff) ? 1.0 : 0.0;
    for (int seq_idx = 0; seq_idx < n_seqs_eval; ++seq_idx) {
      window_prob[seq_idx] = deterministic_window_prob;
    }
    UNPROTECT(1);
    return result;
  }

  // scale factor for the PWM, so the summed integer range stays ~resolution regardless of motif width.
  const double score_bin_scale = resolution / score_range;

  // Convert PWM scores above each position-specific minimum to integer bin shifts.
  // Flooring (round down) keeps discretized scores conservative.
  // 'volatile' forces separate floating-point operations to match R's binning.
  int* pwm_score_shift = (int*) R_alloc(4 * motif_width, sizeof(int)); // score-bin shift for each nt at each motif pos
  int* max_score_shift_by_pos = (int*) R_alloc(motif_width, sizeof(int)); // max score-bin shift at each motif pos
  long long max_reachable_score_bin = 0; // max summed score-bin shift across the full PWM
  for (int mot_pos = 0; mot_pos < motif_width; ++mot_pos) {
    int max_score_shift = 0;
    for (int nt = 0; nt < 4; ++nt) {
      volatile double score_offset = pwm[nt + 4 * mot_pos] - min_score_by_pos[mot_pos];
      volatile double scaled_score_offset = score_offset * score_bin_scale;
      const int score_shift = static_cast<int>(floor_grid(scaled_score_offset));
      pwm_score_shift[nt + 4 * mot_pos] = score_shift;
      if (score_shift > max_score_shift) max_score_shift = score_shift;
    }
    max_score_shift_by_pos[mot_pos] = max_score_shift;
    max_reachable_score_bin += max_score_shift;
  }

  // exclude first score bin strictly above the cutoff through floor(bin) + 1 -> excludes scores equal to the cutoff.
  volatile double cutoff_offset = cutoff - min_pwm_score;
  volatile double scaled_cutoff = cutoff_offset * score_bin_scale;
  const int cutoff_bin = static_cast<int>(floor_grid(scaled_cutoff)) + 1;

  // cutoff outside the reachable score range result in all/no seqs having motif occ.; P(m) = 1 or 0.
  if (cutoff_bin <= 0) {
    for (int seq_idx = 0; seq_idx < n_seqs_eval; ++seq_idx) window_prob[seq_idx] = 1.0;
    UNPROTECT(1);
    return result;
  }
  if (cutoff_bin > max_reachable_score_bin) {
    for (int seq_idx = 0; seq_idx < n_seqs_eval; ++seq_idx) window_prob[seq_idx] = 0.0;
    UNPROTECT(1);
    return result;
  }

  // track only bins reachable after each motif pos, capped at the cutoff bin.
  // the cutoff bin collects all prob mass that has already passed the cutoff.
  int* n_score_bins_by_pos = (int*) R_alloc(motif_width, sizeof(int));
  {
    long long max_reachable_bin = 0;
    for (int mot_pos = 0; mot_pos < motif_width; ++mot_pos) {
      max_reachable_bin += max_score_shift_by_pos[mot_pos];
      const long long max_tracked_bin = std::min<long long>(max_reachable_bin, cutoff_bin);
      n_score_bins_by_pos[mot_pos] = static_cast<int>(max_tracked_bin) + 1;
    }
  }
  const int n_allocated_score_bins = n_score_bins_by_pos[motif_width - 1];

  if (!order1) { // 0th-order: for each seq, track one prob distribution across score bins.
    double* score_prob = (double*) R_alloc(n_allocated_score_bins, sizeof(double));
    double* updated_score_prob = (double*) R_alloc(n_allocated_score_bins, sizeof(double));

    for (int seq_idx = 0; seq_idx < n_seqs_eval; ++seq_idx) {
      if ((seq_idx % INTERRUPT_CHECK_INTERVAL) == 0) R_CheckUserInterrupt();

      int n_score_bins = 1;
      score_prob[0] = 1.0; // before scoring the PWM, all prob mass is in score bin 0

      // construct the summed score distribution one PWM position at a time.
      for (int mot_pos = 0; mot_pos < motif_width; ++mot_pos) {
        const int n_updated_bins = n_score_bins_by_pos[mot_pos];
        std::fill(updated_score_prob, updated_score_prob + n_updated_bins, 0.0);

        // shift the current dist by each nt score and weight it by the nt prob.
        for (int nt = 0; nt < 4; ++nt) {
          const double nt_prob = seqs_freq_nt[seq_idx + n_seqs_eval * nt];
          if (nt_prob == 0.0) continue;
          add_weighted_score_prob(score_prob, n_score_bins, pwm_score_shift[nt + 4 * mot_pos],
                                  nt_prob,updated_score_prob, n_updated_bins);
        }
        // use the updated dist as input at the next PWM position.
        std::swap(score_prob, updated_score_prob);
        n_score_bins = n_updated_bins;
      }
      window_prob[seq_idx] = clamp01(score_prob[n_score_bins - 1]); // last bin holds all passing mass
    }

  } else { // 1st-order: track a separate score distribution for partial windows ending in each nt.
    // track current and updated score distributions for branches ending in A, C, G, or T.
    double* score_prob_by_nt[4];
    double* updated_score_prob_by_nt[4];
    for (int nt = 0; nt < 4; ++nt) {
      score_prob_by_nt[nt] = (double*) R_alloc(n_allocated_score_bins, sizeof(double));
      updated_score_prob_by_nt[nt] = (double*) R_alloc(n_allocated_score_bins, sizeof(double));
    }

    for (int seq_idx = 0; seq_idx < n_seqs_eval; ++seq_idx) {
      if ((seq_idx % INTERRUPT_CHECK_INTERVAL) == 0) R_CheckUserInterrupt();

      // initialize the first PWM pos using the marginal nt probs.
      int n_score_bins = n_score_bins_by_pos[0];
      for (int nt = 0; nt < 4; ++nt) {
        std::fill(score_prob_by_nt[nt], score_prob_by_nt[nt] + n_score_bins, 0.0);
        const int score_bin = std::min(pwm_score_shift[nt + 4 * 0], n_score_bins - 1);
        score_prob_by_nt[nt][score_bin] += seqs_freq_nt[seq_idx + n_seqs_eval * nt];
      }

      // extend each branch using P(current nt | previous nt) at the remaining PWM positions.
      for (int mot_pos = 1; mot_pos < motif_width; ++mot_pos) {
        const int n_updated_bins = n_score_bins_by_pos[mot_pos];
        for (int current_nt = 0; current_nt < 4; ++current_nt) {
          std::fill(
            updated_score_prob_by_nt[current_nt],
                                    updated_score_prob_by_nt[current_nt] + n_updated_bins,
                                    0.0);
          const int nt_score_shift = pwm_score_shift[current_nt + 4 * mot_pos];
          for (int previous_nt = 0; previous_nt < 4; ++previous_nt) {
            const double nt_prob = cond_prob[current_nt + 4 * previous_nt + 16 * seq_idx];
            if (nt_prob == 0.0) continue;
            add_weighted_score_prob(
              score_prob_by_nt[previous_nt], n_score_bins, nt_score_shift, nt_prob,
              updated_score_prob_by_nt[current_nt], n_updated_bins);
          }
        }
        // use the updated distributions as input at the next PWM pos.
        for (int nt = 0; nt < 4; ++nt) {
          std::swap(score_prob_by_nt[nt], updated_score_prob_by_nt[nt]);
        }
        n_score_bins = n_updated_bins;
      }

      // sum the passing prob across all possible ending nucleotides.
      double passing_prob = 0.0;
      for (int nt = 0; nt < 4; ++nt) {
        passing_prob += score_prob_by_nt[nt][n_score_bins - 1];
      }
      window_prob[seq_idx] = clamp01(passing_prob);
    }
  }

  UNPROTECT(1);
  return result;
}
