# Simulation Script
rm(list=ls())
setwd(this.path::here())

library(sparseVARboot)
library(parallel)

#### Sourcing scripts ####
source("DGP.R") # R-script that collects relevant DGPs
source("Simulation.R")
#source("Common_Parameters.R")
################################################################################
# Sample sizes
n <- c(25, 50, 75, 100)
N <- c(20, 40, 100, 200)
# Simulations
sim <- 500
# Bootstrap replications
B <- 199
# Confidence level
level <- c(0.9, 0.95, 0.99)
# Null hypothesis
mu0 <- 0
# Methods
boot <- c(#"VAR-L1-unpen-own-BIC",
          "VAR-L1-unpen-own-TF-11",
          #"MBB", "BWB", "DWB",
          #"VAR-GP-unpen-own-BIC",
          "VAR-GP-unpen-own-TF-11")
# parallel
parallel_sims <- TRUE
# test statistic is absolute value
abs_val <- TRUE
# naive standardization
standardize <- FALSE
################################################################################

set.seed(2003202052)

# DGPs
type <- 2
# Mean
mus <- seq(from = -0.25, to = 0.25, length.out = 51)
# Proportion
prop <- 0.5
# Methods
#boot <- c(boot, "VAR-oracle")

for(m in 1:length(mus)){
  pars <- expand.grid(mean = mus[m], prop = prop, n = n, N = N, DGP = type)
  out <- simulations(pars = pars, boot = boot, mu0 = mu0, sim = sim, B = B, level = level, p = 0, l = 0,
                     abs_val = abs_val, standardize = standardize, parallel_sims = parallel_sims)
  reject <- out$reject
  tuning <- out$tuning
  timing <- out$timing

  save(reject, tuning, timing, file = paste0("../Results/power_curve_results/dgp2_power_curve_mu", m, ".RData"))
}
