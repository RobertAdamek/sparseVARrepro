# simulate_boot <- function(pars, boot, B, level, p = 0, l = 0, selection = 1, parallel_sims = TRUE) {
#   if (parallel_sims) {
#     n_cores <- 1
#   } else {
#     n_cores <- NULL
#   }
#   reject <- array(dim = c(length(boot), length(level)))
#   dimnames(reject) <- list(boot = boot, level = 1 - level)
#   
#   tuning <- rep(NA, length(boot))
#   names(tuning) <- boot
#   
#   x <- sim_DGP(n = pars$n, N = pars$N, type = pars$DGP, mu = pars$mean, prop = pars$prop)
#   for (b in 1:length(boot)) {
#     boot_method <- 1*(boot[b] == "VAR-L1") + 1*(boot[b] == "VAR-HL") + 
#       4*(boot[b] == "MBB") + 3*(boot[b] == "BWB")
#     pen <- 1*(boot[b] == "VAR-L1") + 2*(boot[b] == "VAR-HL")
#     out <- boot_means(x = x, boot = boot_method, penalization = pen, p = p, l = l, 
#                       B = B, q = level, selection = selection, show_progress = FALSE, n_cores = n_cores)
#     reject[b, ] <- out$mean > out$boot_quantiles
#     tuning[b] <- out$par
#   }
#   return(list(reject = reject, tuning = tuning))
# }
# 
# simulate_boot_penalization <- function(pars, boot, B, level, p = 0, l = 0, selection = 1, parallel_sims = TRUE) {
#   if (parallel_sims) {
#     n_cores <- 1
#   } else {
#     n_cores <- NULL
#   }
#   reject <- array(dim = c(length(boot), length(level)))
#   dimnames(reject) <- list(boot = boot, level = 1 - level)
#   
#   tuning <- rep(NA, length(boot))
#   names(tuning) <- boot
#   
#   x <- sim_DGP(n = pars$n, N = pars$N, type = pars$DGP, mu = pars$mean, prop = pars$prop)
#   for (b in 1:length(boot)) {
#     boot_method <- 1 # VAR bootstrap
#     pen <- 1 #VAR bootstrap with L1 penalizqtion 
#     
#     if(boot[b]=="VAR-L1-pen"){
#       pen_own <- TRUE # Penalize own lags
#     }else{
#       pen_own <- FALSE # Do not penalize own lags
#     }
#     if(boot[b]=="VAR-L1-unpen-own-1"){
#       only_lag1 <- TRUE # Only penalize own lag 1
#     }else{
#       only_lag1 <- FALSE # Penalize own lag 1 to p
#     }
# 
#     out <- boot_means(x = x, boot = boot_method, penalization = pen, p = p, l = l, 
#                       B = B, q = level, selection = selection, show_progress = FALSE, n_cores = n_cores,
#                       pen_own = pen_own, only_lag1 = only_lag1)
#     
#     reject[b, ] <- out$mean > out$boot_quantiles
#     tuning[b] <- out$par
#   }
#   return(list(reject = reject, tuning = tuning))
# }

simulate_boot_all_methods <- function(pars, boot, B, level, p = 0, l = 0,
                                      abs_val = FALSE, standardize = FALSE, 
                                      parallel_sims = TRUE) {
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
  selection < 1
  pen <- 1
  pen_own <- TRUE # Penalize own lags
  only_lag1 <- FALSE # Penalize own lag 1 to p
  
  sD <- sim_DGP(n = pars$n, N = pars$N, type = pars$DGP, mu = pars$mean, prop = pars$prop)
  x <- sD$x
  for (b in 1:length(boot)) {
    if (substring(boot[b], 1, 6) == "VAR-L1") {
      boot_method <- 1
      pen <- 1 #VAR bootstrap with L1 penalizqtion 
      if (grepl("BIC", boot[b])) {
        selection <- 1
      } else if (grepl("PI", boot[b])) {
        selection <- 4
        PI_c <- as.numeric(substring(boot[b], nchar(boot[b]) - 2, nchar(boot[b]))) / 10
      } else if (grepl("TF", boot[b])) {
        selection <- 5
        PI_c <- as.numeric(substring(boot[b], nchar(boot[b]) - 2, nchar(boot[b]))) / 10
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
    
    out <- boot_means(x = x, oracle_A=sD$A, oracle_u=sD$u, boot = boot_method, penalization = pen, p = p, l = l, 
                      abs_val, standardize,
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
    
    if (abs_val) {
      reject[b, ] <- abs{out$mean} > out$boot_quantiles
    } else {
      reject[b, ] <- out$mean > out$boot_quantiles
    }
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
