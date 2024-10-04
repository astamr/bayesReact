#' Wrapper function to construct a STAN model for motif activity inference.
#' @description
#' This function constructs a STAN model object for motif activity inference, and is used by fit_motif_model(), bayesReact_core(), and bayesReact_parallel().
#'
#' @param model A character string specifying the model to be used for inference:
#' Either "bayesReact" (default) or "BF" (Bayes Factor; when comparing support for the beta model against the uniform null model).
#'
#' @return A stanmodel object.
#' @export
#'
#' @examples
#' \dontrun{
#' my_model <- construct_motif_model()
#' }
#'
construct_motif_model <- function(model = "bayesReact"){

  ## Model definitions (STAN code) ##
  # Core one-parameter model for motif activity inference
  bayesReact <- 'data {           // default model    ### CONSIDER ADDING M PLATE ###
  int<lower=1> K;                 // number of motif observations
  int<lower=1> C;                 // number of samples or cells

  vector[C] sum_log_l;            // precomputed sum log l for each sample, relevant when adding M plate
  vector[C] sum_log_r;            // precomputed sum log r for each sample
  vector[C] sum_log_1_minus_r;    // precomputed sum log 1-r for each sample
  }
  parameters {
    vector[C] a;                  // core shape parameter for each sample
  }
  model {
    //priors
    a ~ normal(0, 10);

    //likelihood - n|t ~ multinomial
    for (c in 1:C) {
      if (a[c] < 0) {
        target += sum_log_l[c] + (-a[c])*sum_log_1_minus_r[c] - K*lbeta(1, 1-a[c]);
      }
      else {
        target += sum_log_l[c] + a[c]*sum_log_r[c] - K*lbeta(1+a[c], 1);
      };
    }
  }
  '

  # Core one-parameter model for motif activity inference with shrinkage prior (horseshoe)
  bayesReact_shrinkage <- 'data { // shrinkage prior model
  int<lower=1> K;                 // number of motif observations
  int<lower=1> C;                 // number of samples or cells

  vector[C] sum_log_l;            // precomputed sum log l for each sample, relevant when adding M plate
  vector[C] sum_log_r;            // precomputed sum log r for each sample
  vector[C] sum_log_1_minus_r;    // precomputed sum log 1-r for each sample
  }
  parameters {
    real<lower=0> tau;
    vector<lower=0>[C] lambda;
    vector[C] z;
  }
  transformed parameters {
  vector[C] a;                    // core shape parameter for each sample
  real<lower=0> cnst = 0.1;       // Regularization parameter/constant

  a = z .* (tau * lambda) ./ sqrt(1 + square(tau * lambda / cnst));
  }
  model {
    //priors
    tau ~ cauchy(0, 2.5);        // Global shrinkage parameter
    lambda ~ cauchy(0, 2.5);     // Local shrinkage parameter
    z ~ normal(0, 1);            // Unscaled activity parameter (?)

    //likelihood - n|t ~ multinomial
    for (c in 1:C) {
      if (a[c] < 0) {
        target += sum_log_l[c] + (-a[c])*sum_log_1_minus_r[c] - K*lbeta(1, 1-a[c]);
      }
      else {
        target += sum_log_l[c] + a[c]*sum_log_r[c] - K*lbeta(1+a[c], 1);
      };
    }
  }
  '

  # BF: model comparison between beta model and uniform null model
  BF <- 'data {                 // fits model for single sample at a time (C = 1)
  int<lower=1> K;               // number of motif observations

  real sum_log_l;               // precomputed sum log l for each sample, relevant when adding M plate
  real sum_log_r;               // precomputed sum log r for each sample
  real sum_log_1_minus_r;       // precomputed sum log 1-r for each sample
  }
  parameters {
    real a;                     // core shape parameter for each sample
  }
  model {
    //prior
    a ~ normal(0, 10);

    //likelihood - n|t ~ multinomial
    if (a < 0) {
      target += sum_log_l + (-a)*sum_log_1_minus_r - K*lbeta(1, 1-a);
      }
    else {
      target += sum_log_l + a*sum_log_r - K*lbeta(1+a, 1);
    };
  }
  '

  # Core two-parameter model for motif activity inference
  bayesReact_2param <- 'data {           // default model    ### CONSIDER ADDING M PLATE ###
  int<lower=1> K;                 // number of motif observations
  int<lower=1> C;                 // number of samples or cells

  vector[C] sum_log_l;            // precomputed sum log l for each sample, relevant when adding M plate
  vector[C] sum_log_r;            // precomputed sum log r for each sample
  vector[C] sum_log_1_minus_r;    // precomputed sum log 1-r for each sample
  }
  parameters {
    vector<lower=0>[C] alpha;                  // shape1 parameter for each sample
    vector<lower=0>[C] beta;                   // shape2 parameter for each sample
  }
  model {
    //priors
    //a ~ normal(0, 10);
    alpha ~ exponential(0.1); // leads to sd = 10
    beta ~ exponential(0.1);

    //likelihood - n|t ~ multinomial
    for (c in 1:C) {
      target += sum_log_l[c] + (alpha[c]-1)*sum_log_r[c] + (beta[c]-1)*sum_log_1_minus_r[c] - K*lbeta(alpha[c], beta[c]);
      };
    }
  generated quantities {
    vector[C] a;
    for (c in 1:C) {
      a[c] = alpha[c] / (alpha[c] + beta[c]);
    }
  }
  '

  # Extended model accounting for target efficiency (e)
  # ...

  # list of models
  model_codes <- list(bayesReact = bayesReact, BF = BF,
                      bayesReact_2param = bayesReact_2param, bayesReact_shrinkage = bayesReact_shrinkage)

  # Construct and return STAN model
  return(rstan::stan_model(model_code = model_codes[[model]]))
}
