stepdown_robert<-function(boot_means_output, levels=c(0.9,0.95,0.99), hypothesis=0, verbose_output=FALSE){
  bmo<-boot_means_output
  smeans<-bmo$smeans; smeans[,2]<-smeans[,2]+1;smeans<-smeans[order(smeans[,2]), ]
  means_boot<-bmo$means_boot; means_boot[,2,]<-means_boot[,2,]+1
  for(i in 1:dim(means_boot)[3]){
    means_boot[,,i]<-means_boot[order(means_boot[,2,i]),,i]
  }
  means_boot[,1,]<-means_boot[,1,]+hypothesis
  smeans_0<-smeans; means_boot_0<-means_boot
  verbose<-list()
  simple<-list()
  for(ell in 1:length(levels)){
    smeans<-smeans_0; means_boot<-means_boot_0
    verbose[[paste("level",levels[ell])]]=list()
    iterations<-1
    in_set<-1:nrow(bmo$smeans)
    rejected_set<-NULL
    newly_rejected<-1
    while(newly_rejected>0 && length(in_set)>0){
      newly_rejected<-0
      max_boot_means<-apply(means_boot, 3, function(x){max(x[,1])})
      boot_quantile<-quantile(max_boot_means,levels[ell])
      verbose[[paste("level",levels[ell])]][[paste("step",iterations)]]=list("quantile"=boot_quantile, "rejected_means"=NULL,"rejects"=NULL)
      for(i in 1:nrow(smeans)){
        if(smeans[i,1]>boot_quantile){
          verbose[[paste("level",levels[ell])]][[paste("step",iterations)]]$rejected_means<-c(verbose[[paste("level",levels[ell])]][[paste("step",iterations)]]$rejected_means, smeans[i,1]) # add the rejected mean to the output
          rejected_set<-c(rejected_set, smeans[i,2])# update the rejected set and the output
          verbose[[paste("level",levels[ell])]][[paste("step",iterations)]]$rejects<-c(verbose[[paste("level",levels[ell])]][[paste("step",iterations)]]$rejects, smeans[i,2]) # add the indices of the rejected units to the output
          in_set<-in_set[-which(in_set==smeans[i,2])] # update the surviving set
          newly_rejected<-newly_rejected+1 # keep count of how many were rejected
        }
      }
      smeans<-smeans_0[in_set,] # remove the rejected units from the original means
      means_boot<-means_boot_0[in_set,,] # also remove them from the other bootstrapped statistics
      iterations<-iterations+1
    }
    simple[[paste("level", levels[ell], "rejects")]]=rejected_set
  }
  if(verbose_output){
    return(verbose)
  }else{
    return(simple)
  }
}

StepM <- function(x, names_units = NULL, mu0 = 0, boot = 1, p = 0, l = 0,
                  abs_val = FALSE, standardize = TRUE, q = c(0.9, 0.95, 0.99), B = 9999,
                  show_progress = TRUE, penalization = 1,
                  selection = 5, c = 1.1, pen_own = FALSE, only_lag1 = FALSE) {
  # x is the data matrix
  # oracle_A and oracle_u are the true VAR coefficients and errors respectively. Only used when simulating an oracle method, choosing penalization=-1
  # boot indicates which bootstrap method to use. 1=sparse VAR bootstrap, 3=BWB, 4=MBB
  # p is the number of lags used when estimating the VAR. p=0 gives automatic lag selection
  # l is the block length when using one of the block bootstrap methods. l=0 gives automatic block length selection
  # abs_val is a boolean if we should take the absolute value of the test statistic. FALSE for one-sided tests
  # q are the levels of the test
  # B are the bootstrap replications
  # show_progress is a boolean if a progres bar should be shown
  # penalization: integer, 0 (no penalization), 1 (L1), 2 (HLag)
  # nbr_lambdas : double, number of sparsity parameters to consider in grid (for simplicity set as double)
  # lambda_ratio : double, ratio lambda_max/lambda_min
  # selection: integer, 1 for bic, 2 for aic, 3 for hq, 4 for plug-in, 5 for theoretically founded, to select the sparsity parameter lambda
  # eps: double, convergence tolerance
  # pen_own: boolean: true if own lags are to be penalized, false if they should be unpenalized
  # only_lag1 : boolean, only relevant if pen_own = false : TRUE if only first lag should be unpenalized, FALSE all own lags should be unpenalized
  # c is the constant used in the plug-in and theoretically founded selection methods
  # K is the number of iterations in the plug-in and theoretically founded selection methods
  # improvement_thresh: is a threshold to stop iterating the plug-in/TF methods. 0.01 means it stops if there's less than 1% improvement in the objective function
  # Nsim is the number of gaussians drawn when computing the quantiles of the max of correlated gaussians in the plug-in selection method
  # alpha is the probability parameter in the plug-in and TF methods. Lower alpha gives bigger lambdas.
  nN <- ncol(x)
  method_details <- vector(mode = "list", length = length(boot))
  names(method_details) <- boot
  for (b in 1:length(boot)) {
    boot_method <- 1*(boot[b] == "VAR-MB") + 2*(boot[b] == "VAR-EB") + 
      4*(boot[b] == "MBB") + 3*(boot[b] == "BWB") + 3*(boot[b] == "MB") + 4*(boot[b] == "EB")
    if (boot[b] %in% c("MB", "RB")) {
      l <- 1
    }
    out <- boot_means_clean(x = x, mu0 = mu0, boot = boot_method, p = p, l = l,
                          abs_val = abs_val, standardize = standardize, q = q, B = B,
                          show_progress = show_progress, penalization = penalization,
                          selection = selection, c = c, pen_own = pen_own, only_lag1 = only_lag1)
    method_details[[b]] <- out
  
    i <- 1 + out$smeans[, 2]
    t_r <- out$smeans[, 1]
    t_b <- sapply(1:B, function(b){
      tm <- rep(NA, nN)
      tm[1 + out$means_boot[, 2, b]] <- out$means_boot[, 1, b]
      return(tm)
    })
    
    max_t_r_b <- sapply(1:B, function(b){
      rev(cummax(t_b[rev(i), b]))
    })
    p_init <- sapply(1:nN, function(r){
      (sum(max_t_r_b[r, ] >= t_r[r]) + 1) / (B + 1)
    })
    
    p_adj <- cummax(p_init)
    if (b == 1){
      p_val <- data.frame(unit = names_units[i], 
                          mean_growth = colMeans(x)[i])
    }
    p_val[[boot[b]]] <- p_adj
  }
  return(list(p_val = p_val, details = method_details))
}

