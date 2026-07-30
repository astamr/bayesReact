#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <algorithm>
#include <vector>

// This function counts every possible k-mer separately for each sequence, using greedy
// non-overlapping matches (equivalent to Regmex::n.obs.mot(..., overlap = FALSE)).
// Scanning advances one nucleotide at a time. After accepting a k-mer starting
// at position pos, that same k-mer cannot be accepted at starts pos + 1 through
// pos + k - 1; it becomes eligible again at pos + k.

// Only uppercase A/C/G/T in (motifs and) seqs; maps A,C,G,T -> 0,1,2,3.
static inline unsigned int nt_to_code(const unsigned char nt) {
  const unsigned int bits = (nt >> 1) & 3u;
  return bits ^ (bits >> 1);
}

extern "C" SEXP fast_kmers_counting(SEXP seqs, SEXP k_len) {
  const int k = INTEGER(k_len)[0];
  const int n_seqs = XLENGTH(seqs);
  const int n_kmers = 1 << (2 * k);

  // Allocate and protect an n_seqs × n_kmers R matrix, initialise all counts to zero, and obtain direct C++ access to its numeric values.
  SEXP motif_counts = PROTECT(Rf_allocMatrix(REALSXP, n_seqs, n_kmers));
  double* motif_counts_i = REAL(motif_counts);
  std::fill(motif_counts_i,
            motif_counts_i + static_cast<R_xlen_t>(n_seqs) * n_kmers,
            0.0);

  // Block a k-mer from being counted at its next k - 1 overlapping start
  // positions after an accepted match. `release` records when it is eligible
  // to be counted again.
  std::vector<unsigned char> blocked(n_kmers, 0);
  std::vector<int> release(k, -1);
  const unsigned int mask = static_cast<unsigned int>(n_kmers - 1);

  // Scan each sequence independently. Reset the sliding k-mer code, sequence
  // position, and k-mer unblocking index so non-overlap is enforced separately within each sequence.
  for (int seq_index = 0; seq_index < n_seqs; ++seq_index) {
    const char* seq = CHAR(STRING_ELT(seqs, seq_index));
    unsigned int current_kmer_code = 0;
    int seq_pos = 0;
    int kmer_unblock_index = 0;

    // Move through the current sequence one nucleotide at a time.
    for (const unsigned char* p =
         reinterpret_cast<const unsigned char*>(seq);
         *p != '\0'; ++p) {
      // Append the current nt to the sliding kmer-code and discard the nt that has moved outside the current k-mer window.
      current_kmer_code = ((current_kmer_code << 2) | nt_to_code(*p)) & mask;
      if (++seq_pos < k) continue; // A full k-mer cannot be evaluated until k nts are observed.

      // Release the k-mer that was blocked k positions earlier, allowing it to be
      // counted again if it occurs at the current non-overlapping position.
      const int kmer_code_to_release = release[kmer_unblock_index];
      if (kmer_code_to_release >= 0) blocked[kmer_code_to_release] = 0;

      const int kmer_code = static_cast<int>(current_kmer_code); // The current sliding code window represents the full k-mer ending at seq_pos.
      release[kmer_unblock_index] = -1; // Clear the current release record before reuse.

      // Count given k-mer only if it was not counted at an overlapping start position.
      // Then block it until it becomes eligible again.
      if (!blocked[kmer_code]) {
        //motif_counts_i[seq_index + n_seqs * kmer_code] += 1.0; // Update to avoid overflow
        const R_xlen_t matrix_index =
          static_cast<R_xlen_t>(seq_index) +
          static_cast<R_xlen_t>(n_seqs) * kmer_code;
        motif_counts_i[matrix_index] += 1.0;

        blocked[kmer_code] = 1;
        release[kmer_unblock_index] = kmer_code;
      }

      if (++kmer_unblock_index == k) kmer_unblock_index = 0; // Resets release index after k positions.
    }

    // Unblock k-mers still blocked at the end of this sequence before scanning the next sequence.
    for (int unblock_index = 0; unblock_index < k; ++unblock_index) {
      const int kmer_code = release[unblock_index];
      if (kmer_code >= 0) {
        blocked[kmer_code] = 0;
        release[unblock_index] = -1;
      }
    }
  }

  UNPROTECT(1);
  return motif_counts;
}

// Register native functions that R may call via .Call()
static const R_CallMethodDef CallEntries[] = {
  {"fast_kmers_counting",
   (DL_FUNC) &fast_kmers_counting, 2},
   {NULL, NULL, 0}
};

// Called when the bayesReact package DLL is loaded.
extern "C" void R_init_bayesReact(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
