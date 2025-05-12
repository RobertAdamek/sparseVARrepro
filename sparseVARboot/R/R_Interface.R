#' @export
sparseVAR <- function(Y, p, trim, pen, nbr_lambdas, lambda_ratio, eps, selection, c = 0.8, K = 15, improvement_thresh = 0.01, 
                        Nsim = 1000, alpha = 0.05, pen_own = T, only_lag1 = F) {
  sparseVAR_R(Y, p, trim, pen, nbr_lambdas, lambda_ratio, eps, selection, c, K, improvement_thresh, Nsim, alpha, pen_own, only_lag1)
}

#' @export
boot_means <- function(x, oracle_A, oracle_u, boot = 1, p = 0, l = 0, abs_val = TRUE, q = 0.95, B = 9999,
                       penalization = 1, nbr_lambdas = 10, lambda_ratio = 100,
                       selection = 1, eps = 0.001, show_progress = TRUE, n_cores = NULL,
                       pen_own = T, only_lag1 = F, c = 0.8, K = 15, improvement_thresh = 0.01, 
                       Nsim = 1000, alpha = 0.05) {
  if (is.null(n_cores)) {
    n_cores <- parallelly::availableCores(omit = 2)
  }
  RcppParallel::setThreadOptions(numThreads = n_cores)
  
  boot_means_R(x, oracle_A, oracle_u, boot, p, l, abs_val, q, B, show_progress,
               penalization, nbr_lambdas, lambda_ratio, selection, eps, pen_own, only_lag1, c, K, improvement_thresh, Nsim, alpha)
}

#' @export
sim_VAR <- function(n, A, Sigma, burn = 50) {
  sim_VAR_cpp_both(n, A, Sigma, burn)
}