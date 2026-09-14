simulations <- function(pars, boot, mu0 = 0, sim = 1000, B = 199, level = 0.95, p = 0, l = 0, 
                        abs_val = TRUE, standardize = FALSE, parallel_sims = TRUE) {
  parsnames <- paste0("(DGP ", pars$DGP, ", N = ", pars$N, ", n = ", pars$n, 
                      ", prop = ", pars$prop, ", mu = ", pars$mean, 
                      ")")
  
  reject <- array(dim = c(nrow(pars), length(boot), length(level)))
  dimnames(reject) <- list(pars = parsnames, boot = boot, level = 1 - level)
  
  tuning <- array(dim = c(nrow(pars), length(boot), sim))
  dimnames(tuning) <- list(pars = parsnames, boot = boot, sim = 1:sim)
  
  seeds <- sample.int(2^20, size = sim)
  
  if (parallel_sims) {
    cl <- parallel::makeCluster(parallelly::availableCores(omit = 2))
    parallel::clusterExport(cl, varlist = ls(globalenv()))
    parallel::clusterEvalQ(cl, library(sparseVARboot))
    parallel::clusterSetRNGStream(cl, sample.int(2^20, size = 1))
  }
  
  for (i in 1:nrow(pars)) {
    parsi = pars[i, ]
    if (parallel_sims) {
      out <- parallel::parLapply(cl = cl, X = seeds, 
                                 fun = simulate_boot_all_methods, 
                                 pars = pars[i, ], mu0 = mu0, boot = boot, B = B, level = level, p = 0, l = 0, 
                                 abs_val = abs_val, standardize = standardize, parallel_sims = parallel_sims)
    } else {
      out <- lapply(X = seeds, 
                    FUN = simulate_boot_all_methods, 
                    pars = pars[i, ], mu0 = mu0, boot = boot, B = B, level = level, p = 0, l = 0, 
                    abs_val = abs_val, standardize = standardize, parallel_sims = parallel_sims)
    }
    reject[i, , ] <- apply(sapply(out, function(x){x$reject}, simplify = "array"), 1:2, mean)
    tuning[i, , ] <- sapply(out, function(x){x$tuning}, simplify = "array")
  }
  
  if (parallel_sims) {
    parallel::stopCluster(cl)
  }
  return(reject)
}

simulate_boot_all_methods <- function(seed, pars, mu0, boot, B, level, p = 0, l = 0,
                                      abs_val = TRUE, standardize = FALSE, 
                                      parallel_sims = TRUE) {
  set.seed(seed)
  if (parallel_sims) {
    n_cores <- 1
  } else {
    n_cores <- NULL
  }
  reject <- array(dim = c(length(boot), length(level)))
  dimnames(reject) <- list(boot = boot, level = 1 - level)
  
  tuning <- rep(NA, length(boot))
  ####################remove
  lambdas <- rep(NA, length(boot))
  lambdass <- list()
  coef_pre <- list()
  coef_post <- list()
  ##########################
  names(tuning) <- boot
  
  PI_c <- 1
  selection <- 1
  pen <- 1
  pen_own <- TRUE # Penalize own lags
  only_lag1 <- FALSE # Penalize own lag 1 to p
  
  sD <- sim_DGP(n = pars$n, N = pars$N, type = pars$DGP, mu = pars$mean, prop = pars$prop)
  x <- sD$x
  for (b in 1:length(boot)) {
    if (substring(boot[b], 1, 6) == "VAR-L1" | substring(boot[b], 1, 6) == "VAR-GP") {
      if (substring(boot[b], 1, 6) == "VAR-L1") {
        boot_method <- 1
      } else if (substring(boot[b], 1, 6) == "VAR-GP") {
        boot_method <- 6
      }
      pen <- 1 #VAR bootstrap with L1 penalizqtion 
      if (grepl("BIC", boot[b])) {
        selection <- 1
      } else if (grepl("PI", boot[b])) {
        selection <- 4
        PI_c <- as.numeric(substring(boot[b], nchar(boot[b]) - 1, nchar(boot[b]))) / 10
      } else if (grepl("TF", boot[b])) {
        selection <- 5
        PI_c <- as.numeric(substring(boot[b], nchar(boot[b]) - 1, nchar(boot[b]))) / 10
      }
      if (grepl("-pen-", boot[b])) {
        pen_own <- TRUE # Penalize own lags
      } else {
        pen_own <- FALSE # Do not penalize own lags
      }
      if (grepl("-unpen-own-1", boot[b])) {
        only_lag1 <- TRUE # Only penalize own lag 1
      } else {
        only_lag1 <- FALSE # Penalize own lag 1 to p
      }
    } else if (boot[b] == "VAR-oracle") {
      boot_method <- 1
      pen <- -1 # -1 skips estimation completely
    } else {
      boot_method <- 3*(boot[b] == "BWB") + 4*(boot[b] == "MBB") + 5*(boot[b] == "DWB")
    }
    
    out <- boot_means(x = x, oracle_A = sD$A, oracle_u = sD$u, mu0 = mu0, boot = boot_method, penalization = pen, p = p, l = l, 
                      abs_val = abs_val, standardize = standardize,
                      B = B, q = level, selection = selection, show_progress = FALSE, n_cores = n_cores,
                      pen_own = pen_own, only_lag1 = only_lag1, c = PI_c, 
                      K = 15, improvement_thresh = 0.01, Nsim = 1000, alpha = 0.05)
    
    ####################remove
    coef_pre[[b]] <- out$coef_pre
    coef_post[[b]] <- out$coef_post
    if (length(out$lambda) > 0){
      lambdas[b] <- out$lambda
    }
    if (length(out$lambdas) > 0){
      lambdass[[b]] <- out$lambdas
    }
    ##########################
    
    reject[b, ] <- out$mean > out$boot_quantiles
    tuning[b] <- out$par
  }
  return(list(reject = reject, tuning = tuning
              ###################remove
              , details = list(
                  coef_pre = coef_pre, 
                  coef_post = coef_post,
                  lambda = lambdas, 
                  lambdas = lambdass, 
                  boot_quantiles = out$boot_quantiles, 
                  statistic = out$mean, 
                  smeans = out$smeans,
                  x = x)
              ###########################
              ))
}
