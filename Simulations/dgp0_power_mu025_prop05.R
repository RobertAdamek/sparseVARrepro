# Simulation Script
rm(list=ls())
setwd(this.path::here())

library(sparseVARboot)
library(parallel)

#### Sourcing scripts ####
source("DGP.R") # R-script that collects relevant DGPs
source("Simulation.R")
source("Common_Parameters.R")

set.seed(200320250)

# DGPs
type <- 0 
# Mean
mu <- 0.25
# Proportion
prop <- 0.5

pars <- expand.grid(mean = mu, prop = prop, n = n, N = N, DGP = type)
parsnames <- paste0("(DGP ", pars$DGP, ", N = ", pars$N, ", n = ", pars$n, 
                    ", N_mu = ", ceiling(pars$prop * pars$N), ", mu = ", pars$mean, ")")

reject <- array(dim = c(nrow(pars), length(boot), length(level)))
dimnames(reject) <- list(pars = parsnames, boot = boot, level = 1 - level)

tuning <- array(dim = c(nrow(pars), length(boot), sim))
dimnames(tuning) <- list(pars = parsnames, boot = boot, sim = 1:sim)

seeds<-sample.int(2^20, size = sim)

parallel_sims <- TRUE

if (parallel_sims) {
  cl <- parallel::makeCluster(parallelly::availableCores(omit = 2))
  parallel::clusterExport(cl, varlist = ls(globalenv()))
  parallel::clusterEvalQ(cl, library(sparseVARboot))
  parallel::clusterEvalQ(cl, source("/home/DGP.R"))
  parallel::clusterEvalQ(cl, source("/home/Simulation.R"))
  parallel::clusterSetRNGStream(cl, sample.int(2^20, size = 1))
}

for (i in 1:nrow(pars)) {
  parallel::clusterExport(cl, varlist = "i")
  if (parallel_sims) {
    out <- parallel::parLapply(cl=cl, X=seeds, 
                               fun=simulate_boot_all_methods, 
                               pars=pars[i, ], boot=boot, B=B, level=level, p = 0, l = 0, parallel_sims = parallel_sims
    )
  } else {
    out <- lapply(X=seeds, 
                  FUN=simulate_boot_all_methods, 
                  pars=pars[i, ], boot=boot, B=B, level=level, p = 0, l = 0, parallel_sims = parallel_sims
    )
  }
  reject[i, , ] <- apply(sapply(out, function(x){x$reject}, simplify = "array"), 1:2, mean)
  tuning[i, , ] <- sapply(out, function(x){x$tuning}, simplify = "array")
}

if (parallel_sims) {
  parallel::stopCluster(cl)
}

save(reject, file = "../Results/dgp0_power_mu025_prop05.RData")