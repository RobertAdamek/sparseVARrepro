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

simulate_boot_all_methods <- function(seed ,pars, boot, B, level, p = 0, l = 0, parallel_sims = TRUE) {
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
  
  sD <- sim_DGP(n = pars$n, N = pars$N, type = pars$DGP, mu = pars$mean, prop = pars$prop)
  x <- sD$x
  for (b in 1:length(boot)) {
    if(boot[b] == "VAR-oracle"){
      pen <- -1 # -1 skips estimation completely
    }else{
      pen <- 1 #VAR bootstrap with L1 penalizqtion 
    }
    
    
    boot_method <- 1*(boot[b] == "VAR-L1-pen-BIC") + 1*(boot[b] == "VAR-L1-unpen-own-1-BIC") + 1*(boot[b] == "VAR-L1-unpen-own-BIC") +
                    1*(boot[b] == "VAR-L1-pen-PI-04") + 1*(boot[b] == "VAR-L1-unpen-own-1-PI-04") + 1*(boot[b] == "VAR-L1-unpen-own-PI-04") + 
                      1*(boot[b] == "VAR-L1-pen-PI-08") + 1*(boot[b] == "VAR-L1-unpen-own-1-PI-08") + 1*(boot[b] == "VAR-L1-unpen-own-PI-08") +
                        1*(boot[b] == "VAR-L1-pen-TF-11") + 1*(boot[b] == "VAR-L1-unpen-own-1-TF-11") + 1*(boot[b] == "VAR-L1-unpen-own-TF-11") +
                          1*(boot[b] == "VAR-L1-pen-TF-08") + 1*(boot[b] == "VAR-L1-unpen-own-1-TF-08") + 1*(boot[b] == "VAR-L1-unpen-own-TF-08") +
                            1*(boot[b] == "VAR-L1-pen-TF-04") + 1*(boot[b] == "VAR-L1-unpen-own-1-TF-04") + 1*(boot[b] == "VAR-L1-unpen-own-TF-04") +
                        4*(boot[b] == "MBB") + 3*(boot[b] == "BWB") + 1*(boot[b] == "VAR-oracle")
    selection <- 1*(boot[b] == "VAR-L1-pen-BIC") + 1*(boot[b] == "VAR-L1-unpen-own-1-BIC") + 1*(boot[b] == "VAR-L1-unpen-own-BIC") +
                    4*(boot[b] == "VAR-L1-pen-PI-04") + 4*(boot[b] == "VAR-L1-unpen-own-1-PI-04") + 4*(boot[b] == "VAR-L1-unpen-own-PI-04") +
                      4*(boot[b] == "VAR-L1-pen-PI-08") + 4*(boot[b] == "VAR-L1-unpen-own-1-PI-08") + 4*(boot[b] == "VAR-L1-unpen-own-PI-08") +
                        5*(boot[b] == "VAR-L1-pen-TF-11") + 5*(boot[b] == "VAR-L1-unpen-own-1-TF-11") + 5*(boot[b] == "VAR-L1-unpen-own-TF-11") + 
                          5*(boot[b] == "VAR-L1-pen-TF-08") + 5*(boot[b] == "VAR-L1-unpen-own-1-TF-08") + 5*(boot[b] == "VAR-L1-unpen-own-TF-08") +
                            5*(boot[b] == "VAR-L1-pen-TF-04") + 5*(boot[b] == "VAR-L1-unpen-own-1-TF-04") + 5*(boot[b] == "VAR-L1-unpen-own-TF-04")
    PI_c <- 0.4*(boot[b] == "VAR-L1-pen-PI-04") + 0.4*(boot[b] == "VAR-L1-unpen-own-1-PI-04") + 0.4*(boot[b] == "VAR-L1-unpen-own-PI-04") +
              0.8*(boot[b] == "VAR-L1-pen-PI-08") + 0.8*(boot[b] == "VAR-L1-unpen-own-1-PI-08") + 0.8*(boot[b] == "VAR-L1-unpen-own-PI-08") +
                1.1*(boot[b] == "VAR-L1-pen-TF-11") + 1.1*(boot[b] == "VAR-L1-unpen-own-1-TF-11") + 1.1*(boot[b] == "VAR-L1-unpen-own-TF-11") +
                  0.8*(boot[b] == "VAR-L1-pen-TF-08") + 0.8*(boot[b] == "VAR-L1-unpen-own-1-TF-08") + 0.8*(boot[b] == "VAR-L1-unpen-own-TF-08") +
                    0.4*(boot[b] == "VAR-L1-pen-TF-04") + 0.4*(boot[b] == "VAR-L1-unpen-own-1-TF-04") + 0.4*(boot[b] == "VAR-L1-unpen-own-TF-04") 
    
    if(boot[b]=="VAR-L1-pen-BIC" | boot[b]=="VAR-L1-pen-PI-04" | boot[b]=="VAR-L1-pen-PI-08" | boot[b] == "VAR-L1-pen-TF-11" | boot[b] == "VAR-L1-pen-TF-08" | boot[b] == "VAR-L1-pen-TF-04"){
      pen_own <- TRUE # Penalize own lags
    }else{
      pen_own <- FALSE # Do not penalize own lags
    }
    
    if(boot[b]=="VAR-L1-unpen-own-1-BIC" | boot[b]=="VAR-L1-unpen-own-1-PI-04" | boot[b]=="VAR-L1-unpen-own-1-PI-08" | boot[b] == "VAR-L1-unpen-own-1-TF-11" | boot[b] == "VAR-L1-unpen-own-1-TF-08" | boot[b] == "VAR-L1-unpen-own-1-TF-04"){
      only_lag1 <- TRUE # Only penalize own lag 1
    }else{
      only_lag1 <- FALSE # Penalize own lag 1 to p
    }
    
    
    out <- boot_means(x = x, oracle_A=sD$A, oracle_u=sD$u, boot = boot_method, penalization = pen, p = p, l = l, 
                      B = B, q = level, selection = selection, show_progress = FALSE, n_cores = n_cores,
                      pen_own = pen_own, only_lag1 = only_lag1,  c = PI_c, 
                      K = 15, improvement_thresh = 0.01, Nsim = 1000, alpha = 0.05)
    
    ####################remove
    coef_pre[[b]] <- out$coef_pre
    coef_post[[b]] <- out$coef_post
    if(length(out$lambda)>0){
      lambdas[b] <- out$lambda
    }
    if(length(out$lambdas)>0){
      lambdass[[b]] <- out$lambdas
    }
    ##########################
    
    reject[b, ] <- out$mean > out$boot_quantiles
    tuning[b] <- out$par
  }
  return(list(reject = reject, tuning = tuning
              ###################remove
              , coef_pre=coef_pre, coef_post=coef_post,
              lambda=lambdas, lambdas=lambdass, boot_quantiles=out$boot_quantiles, statistic=out$mean, smeans=out$smeans
              ###########################
              ))
}
