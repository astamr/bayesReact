#' Fit model for motif activity prediction
#'
#' @description Function for predicting motif activity from ranked sequence data using a simple Bayesian model implemented in STAN.
#'
#' @param input path to input file (.rds format) or list object produced by bayesReact::prep_model_input().
#' @param model STAN model object to be used for motif activity inference. The stanmodel object is constructed by construct_motif_model().
#' @param model_type a character string specifying the type of stanmodel object constructed by construct_motif_model():
#' Either "bayesReact" (default) or "BF" (Bayes Factor; when comparing beta model against uniform null model).
#' @param output_type type of output to be returned; either "activity" (only outputs activity score); "activity_summary" (default; outputs activity score, posterior mean and sd for the underlying activity parameter 'a', log P(|a| <= 0 | data), credible intervals, and simple model diagnostics);
#' "full_posterior" (outputs all MCMC iterations after warm-up period); or "full_model" (returns full stanfit model object).
#' @param CI credible interval (CI) to be returned for activity estimates. Default is 80% CI: c(0.10, 0.90).
#' @param iterations number of iterations to be run by the MCMC sampler.
#' @param chains number of independent MCMC chains to be run.
#' @param warmup initial iterations to be discarded for each chain as warm-up/burn-in.
#' @param cores number of cores, default is equal to the number of chains (which is the maximum number of cores that can be utilized by STAN's MCMC sampler). Alternatively, consider parallel::detectCores().
#' @param keep_warmup whether to keep the warm-up iterations or not.
#'
#' @return Activity estimates in a format specified by output_type.
#' @export
#'
#' @examples
#' \dontrun{
#' my_activity <- fit_motif_model(my_data)
#' # returns df with activity estimates, CIs, and simple diagnostics
#' }
#'
fit_motif_model <- function(input, model, model_type = "bayesReact", output_type = "activity_summary", CI = c(0.10, 0.90), #99% CI = c(0.005, 0.995),
                            iterations = 3000, chains = 3, warmup = 500, cores = chains, keep_warmup = F){
  # Try solving Segfault issue
  #try(reg.finalizer(m), silent = T)
  #if (model_type == "BF"){invisible(gc(reset = TRUE, full = TRUE))}

  ## Load data if file_path is provided ##
  if (is.character(input)) {
    input <- readRDS(input)
  }

  ## Model fit (MCMC sampling) ##
  pars = c("a")

  m <- rstan::sampling(model,
                       data = input,
                       iter = iterations,
                       warmup = warmup,
                       chains = chains,
                       cores = cores,
                       pars = pars,
                       save_warmup = keep_warmup,
                       #seed = 42,
                       control = list(adapt_delta = 0.9, # 0.95 (consider, slower)
                                      max_treedepth = 10), # 12 (consider)
                       open_progress = FALSE,
                       show_messages = FALSE,
                       verbose = FALSE, refresh = 0)

  # Handle Segfault issues generated from STAN's memory allocation problem by re-running the MCMC sampler
  times <- 1
  while (is.null(m)) {
    if (times > 10) {
      stop("MCMC sampling failed 10 times. Consider setting MCMC_cores = 1.")
    }
    print("MCMC sampling failed. Trying again.")
    invisible(gc(reset = TRUE, full = TRUE))
    Sys.sleep(1)
    m <- rstan::sampling(model,
                         data = input,
                         iter = iterations,
                         warmup = warmup,
                         chains = chains,
                         cores = cores,
                         pars = pars,
                         save_warmup = keep_warmup,
                         #seed = 42,
                         control = list(adapt_delta = 0.9, # 0.95 (slower)
                                        max_treedepth = 10), # 12 (slower)
                         open_progress = FALSE,
                         show_messages = FALSE,
                         verbose = FALSE, refresh = 0)
    times <- times + 1
  }

  if (model_type != "BF"){cat("Succesfull model fit. \n")}

  ## If "BF", use bridge sampling to obtain log-marginal likelihoods ##
  if (model_type == "BF") {
    m_logml <- bridgesampling::bridge_sampler(m, silent = T, use_neff = F)$logml
    uni_logml <- input$sum_log_l #- input$K*lbeta(1, 1) # Note, there is no model parameters to marginalize and the beta(1, 1) density is utilized.
    lBF <- m_logml - uni_logml # log BF
  }

  ## Return model fit ##
  if (output_type == "activity_summary" | output_type == "activity") { # default
    m_fit_stats <- rstan::summary(m, probs = c(CI[1], CI[2]))$summary
    m_fit_stats <- as.data.frame(m_fit_stats)
    m_fit_stats <- m_fit_stats[-dim(m_fit_stats)[[1]],] # remove last row with lp__
    m_fit_stats$post_prob <- -(stats::pnorm(0, mean = abs(m_fit_stats$mean), sd = m_fit_stats$sd, log.p = T) + log(2))
    m_fit_stats$activity <- m_fit_stats$post_prob*sign(m_fit_stats$mean)

    if (model_type == "BF"){
      m_fit_stats$lBF <- lBF
      return(m_fit_stats[,c("mean", "activity", "post_prob", "sd", paste0(CI[1]*100, "%"), paste0(CI[2]*100, "%"), "n_eff", "Rhat", "lBF")]) # returns data-frame with mean posterior activity estimates, log BF, etc.
    }
    if (output_type == "activity"){
      return(m_fit_stats$activity)
    }
    return(m_fit_stats[,c("mean", "activity", "post_prob", "sd", paste0(CI[1]*100, "%"), paste0(CI[2]*100, "%"), "n_eff", "Rhat")]) # returns data-frame with mean posterior activity estimates and credible intervals

  } else if (output_type == "full_posterior") { # Only for single motif
    return(rstan::extract(m, permuted = TRUE, inc_warmup = FALSE, include = TRUE)) # returns list object with MCMC iterations for each parameter and log-likelihood (lp__)

    } else if (output_type == "full_model") { # Only for single motif
    return(m) # returns S4 class stanfit object

  } else {
    stop("Invalid output_type argument. Please choose from 'activity', 'activity_summary', 'full_posterior', or 'full_model'.")
  }
}







