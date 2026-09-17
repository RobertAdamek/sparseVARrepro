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
mu <- 0
# Proportion
prop <- 1

pars <- expand.grid(mean = mu, prop = prop, n = n, N = N, DGP = type)
out <- simulations(pars = pars, boot = boot, mu0 = mu0, sim = sim, B = B, level = level, p = 0, l = 0, 
                   abs_val = abs_val, standardize = standardize, parallel_sims = parallel_sims)
reject <- out$reject
tuning <- out$tuning
timing <- out$timing

save(reject, tuning, timing, file = "../Results/dgp0_size.RData")