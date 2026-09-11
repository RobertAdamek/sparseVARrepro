# Sample sizes
n <- c(50, 100) # c(50, 100, 200, 500) 
N <- c(20, 40) # c(20, 40, 100, 200)
# Simulations
sim <- 1000
# Bootstrap replications
B <- 199
# Confidence level
level <- c(0.9, 0.95, 0.99)
# Null hypothesis
mu0 <- 0
# Methods
boot <- c("VAR-L1-unpen-own-BIC", 
          "VAR-L1-unpen-own-TF-11", 
          "MBB", "BWB", "DWB",
          "VAR-GP-unpen-own-BIC",
          "VAR-GP-unpen-own-TF-11")
# parallel
parallel_sims <- TRUE
# test statistic is absolute value
abs_val = TRUE
# naive standardization
standardize = FALSE