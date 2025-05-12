#### DGP VAR(1) Krampe et al. (2021), see Appendix D ####
#### https://arxiv.org/pdf/1806.11083.pdf            ####
DGP_Krampe <- function(N, xi = 0.6){
  #### Function to generate data from a sparse VAR(1) process according to Krampe et al. (2021) 
  # N: Number of variables should be either 20, 100 or 200! N needs to be a multiple of 20 for code to work
  # xi: Should be either 0.6, 0.9
  
  # Size of cluster
  block_N <- 20
  
  #### AR parameter matrix ####
  D <- diag(c(xi, -0.7, xi, -0.6, 0.6, 0, 0, 0, 0, 0.2, 0.5, -0.8, 0, 0))
  B <- matrix(c(0.8, 0.2, -0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.6, -0.7, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, -0.9, 0.0, 0.0, 0.0, 0.0, 0.0, -0.6, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.7, 0.0, 0.0, 0.0, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.7,
                0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0),
              nrow = 6, ncol = 14, byrow = T)
  C <- matrix(c(xi, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.3, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, -0.3, 0.0,
                0.6, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.6, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, xi), nrow = 6, ncol = 6, byrow = T)
  A_xi_block <- rbind(cbind(D, matrix(0, nrow = nrow(D), ncol = ncol(C))),
                      cbind(B, C))
  # View(A_xi_block)
  A <- kronecker(diag(1, N/block_N), A_xi_block)
  # View(A)
  
  max_eigenvalue <- max(abs(eigen(A)$values)) 
    # max is 0.8 for xi=0.6, their paper mentions 0.7?
    # max is 0.9 for xi=0.9
  
  
  ### Error Variance Covariance Matrix ####
  Sigma_22 <- matrix(0, 6, 6)
  Sigma_22[, 1] <- Sigma_22[1, ] <- 0.25
  diag(Sigma_22) <- 1
  
  Sigma_11  <- matrix(0, 14, 14)
  diag(Sigma_11[1:5, 1:5][-1, ]) = 0.5
  diag(Sigma_11[1:5, 1:5][, -1]) = 0.5
  diag(Sigma_11[10:12, 10:12][-1, ]) = -0.5
  diag(Sigma_11[10:12, 10:12][, -1]) = -0.5
  diag(Sigma_11) <- 1
  
  Sigma_block <- rbind(cbind(Sigma_11, matrix(0, nrow = nrow(Sigma_11), ncol = ncol(Sigma_22))),
                       cbind(matrix(0, nrow= nrow(Sigma_22), ncol = ncol(Sigma_11)), Sigma_22))
  # View(Sigma_block)
  
  Sigma <- kronecker(diag(1, N/block_N), Sigma_block)
  
  
  #### Output ####
  # A: VAR(1) parameter matrix
  # Sigma: variance-covariance matrix
  # lambda_max: max eigenvalue of A
  
  out <- list("A" = A, "Sigma" = Sigma, "lambda_max" = max_eigenvalue)
}


#### DGPs Kock and Callot (2015)                       ####
#### https://doi.org/10.1016/j.jeconom.2015.02.013     ####

DGP_Kock_A <- function(N){
  #### Function to generate data from a sparse VAR(1) process according to Kock and Callot (2015), Experiment A
  # N: Number of variables should be either 10, 20, 50. 
  
  Phi1 <- diag(0.5, N)
  max_eigenvalue <- max(abs(eigen(Phi1)$values)) 
  Sigma <- diag(.01, N)
  
  #### Output ####
  # Phi1: VAR(1) parameter matrices
  # Sigma: variance-covariance matrix
  # lambda_max: max eigenvalue of companion matrix
  
  out <- list("Phi1" = Phi1,"Sigma" = Sigma, "lambda_max" = max_eigenvalue)
}

DGP_Kock_B <- function(N){
  #### Function to generate data from a sparse VAR(4) process according to Kock and Callot (2015), Experiment B
  # We used this in Adamek, Smeekes, Wilms (2022, JOE)
  # N: Number of variables should be either 10, 20, 50. ! N needs to be a multiple of 5 for code to work
  
  Phi1 <- kronecker(diag(1, N/5), matrix(0.15, 5, 5))
  Phi4 <- kronecker(diag(1, N/5), matrix(-0.1, 5, 5))
  Phi2 <- Phi3 <- matrix(0, N, N)
  companion <- rbind(cbind(Phi1, Phi2, Phi3, Phi4), cbind(diag(1, N*3), matrix(0, 3*N, N)))
  max_eigenvalue <- max(abs(eigen(companion)$values)) 
  Sigma <- diag(.01, N)
  
  #### Output ####
  # Phi1 ... Phi4: VAR(4) parameter matrices
  # companion: matrix companion form
  # Sigma: variance-covariance matrix
  # lambda_max: max eigenvalue of companion matrix
  
  out <- list("Phi1" = Phi1, "Phi2" = Phi2, "Phi3" = Phi3, "Phi4" = Phi4, 
              "companion" = companion, "Sigma" = Sigma, "lambda_max" = max_eigenvalue)
}

DGP_Kock_C <- function(N){
  #### Function to generate data from a sparse VAR(5) process according to Kock and Callot (2015), Experiment C
  # N: Number of variables should be either 10, 20, 50.
  
  Phi1 <- diag(0.95, N)
  Phi2 <- ((-0.95)^(2-1))*Phi1
  Phi3 <- ((-0.95)^(3-1))*Phi1
  Phi4 <- ((-0.95)^(4-1))*Phi1
  Phi5 <- ((-0.95)^(5-1))*Phi1
  
  companion <- rbind(cbind(Phi1, Phi2, Phi3, Phi4, Phi5), cbind(diag(1, N*4), matrix(0, 4*N, N)))
  max_eigenvalue <- max(abs(eigen(companion)$values)) 
  # Is 0.95 but in their paper 0.92...
  
  Sigma <- diag(.01, N)
  
  #### Output ####
  # Phi1 ... Phi5: VAR(5) parameter matrices
  # companion: matrix companion form
  # Sigma: variance-covariance matrix
  # lambda_max: max eigenvalue of companion matrix
  
  out <- list("Phi1" = Phi1, "Phi2" = Phi2, "Phi3" = Phi3, "Phi4" = Phi4, "Phi5" = Phi5, 
              "companion" = companion, "Sigma" = Sigma, "lambda_max" = max_eigenvalue)
}

DGP_Kock_D <- function(N, rho = 0.4){
  #### Function to generate data from a weakly sparse VAR(1) process according to Kock and Callot (2015), Experiment D
  # N: Number of variables should be either 10, 20, 50.
  # rho: 0.4, diag elements of Phi1
  
  Phi1 <- matrix(0, N, N)
  for(i in 1:N){
    for(j in 1:N){
      Phi1[i, j] = ((-1)^(abs(i-j)))*(rho^(abs(i-j) + 1))
    }
  }
  max_eigenvalue <- max(abs(eigen(Phi1)$values)) 
  Sigma <- diag(.01, N)
  
  #### Output ####
  # Phi1: VAR(1) parameter matrix
  # Sigma: variance-covariance matrix
  # lambda_max: max eigenvalue of companion matrix
  
  out <- list("Phi1" = Phi1, "Sigma" = Sigma, "lambda_max" = max_eigenvalue)
}

DGP_single_offdiagonal<-function(N, diagonal=0.8, offdiagonal=0.2){
  Phi1 <- diagonal*diag(N)
  for(i in 2:N){
    Phi1[i-1,i]<-offdiagonal
  }
  max_eigenvalue <- max(abs(eigen(Phi1)$values)) 
  Sigma <- diag(.01, N)
  out <- list("Phi1" = Phi1, "Sigma" = Sigma, "lambda_max" = max_eigenvalue)
}

DGP_special<-function(N, diagonal=0.5, offdiagonal=0.5){
  Phi1 <- diagonal*diag(N)
  for(i in 2:N){
    Phi1[i-1,i]<-offdiagonal
  }
  max_eigenvalue <- max(abs(eigen(Phi1)$values)) 
  Sigma <- toeplitz(c(1,-0.5, rep(0,N-2)))
  out <- list("Phi1" = Phi1, "Sigma" = Sigma, "lambda_max" = max_eigenvalue)
}

sim_DGP_VAR <- function(n, N, type = 1, nburn = 50, ...) {
  if (type == 1) {
    varpars <- DGP_Krampe(N = N, ...)
    A <- t(varpars$A)
    Sigma <- varpars$Sigma
  } else if (type == 2) {
    varpars <- DGP_Kock_A(N = N, ...)
    A <- t(varpars$Phi1)
    Sigma <- varpars$Sigma
  } else if (type == 3) {
    varpars <- DGP_Kock_B(N = N, ...)
    A <- t(cbind(varpars$Phi1, varpars$Phi2, varpars$Phi3, varpars$Phi4))
    Sigma <- varpars$Sigma
  } else if (type == 4) {
    varpars <- DGP_Kock_C(N = N, ...)
    A <- t(cbind(varpars$Phi1, varpars$Phi2, varpars$Phi3, varpars$Phi4))
    Sigma <- varpars$Sigma
  } else if (type == 5) {
    varpars <- DGP_Kock_D(N = N, ...)
    A <- t(varpars$Phi1)
    Sigma <- varpars$Sigma
  } else if (type == 6) {
    varpars <- DGP_Kock_D(N = N, rho=0.2) #Kock D with smaller elements - maximum eigenvalue just below 0.3, so much less persistent than the 0.93 of the default
    A <- t(varpars$Phi1)
    Sigma <- varpars$Sigma
  } else if (type == 7) {
    varpars <- DGP_single_offdiagonal(N = N, ...) #messed up DGP with values on a single off-diagonal
    A <- t(varpars$Phi1)
    Sigma <- varpars$Sigma
  } else if (type == 8) {
    varpars <- DGP_single_offdiagonal(N = N, diagonal=0.3, offdiagonal = 0.7,...) #messed up DGP with values on a single off-diagonal
    A <- t(varpars$Phi1)
    Sigma <- varpars$Sigma
  } else if (type == 9) {
    varpars <- DGP_Kock_D(N = N, rho=0.3) #Kock D with smaller elements - something between DGP 5 and 6
    A <- t(varpars$Phi1)
    Sigma <- varpars$Sigma
  } else if (type == 10) {
    varpars <- DGP_special(N = N,...) #messed up DGP with values on a single off-diagonal. also with special error covariance
    A <- t(varpars$Phi1)
    Sigma <- varpars$Sigma
  }
  
  sV<- sim_VAR(n, A, Sigma, nburn)
  x <- sV$x
  u <- sV$u
  return(list("x" = x, "A" = A, "Sigma" = Sigma, "u"=u))
}

# DGP from Barigozzi, Cho and Owens (2022), "FNETS: Factor-adjusted network estimation 
# and forecasting for high-dimensional time series",
# Appendix E, setting (E.1) - (C.2)
sim_DGP_FNETS <- function(n, N, a = 0.275, 
                          q = 2, burn = 50) {
  U <- matrix(runif(N^2), nrow = N)
  A <- a * (U < 1/N)
  max_ev_A <- max(abs(eigen(A)$values))
  if (max_ev_A > 0.9) {
    A <- 0.9 * A / max_ev_A
  }
  
  D0_diag <- runif(q, min = 0.5, max = 0.8)
  D0 <- matrix(runif(q^2, min = 0, max = 0.3), nrow = q)
  diag(D0) <- D0_diag
  max_ev_D <- max(abs(eigen(D0)$values))
  D <- 0.7 * D0 / max_ev_D
  Lambda <- array(rnorm(2 * N * q), dim = c(N, q, 2))
  
  u = matrix(rnorm((n + burn) * q), ncol = q)
  f = matrix(0, nrow = n + burn, ncol = q)
  for (t in 2:(n + burn)) {
    f[t, ] = f[t-1, ] %*% t(D) + u[t, ]
  }
  chi <- f[-1, ] %*% t(Lambda[, , 1]) + f[-(n + burn), ] %*% t(Lambda[, , 2])
  chi <- chi[1:n + burn - 1, ]
  
  e = matrix(rnorm((n + burn) * N), ncol = N)
  xi = matrix(0, nrow = n + burn, ncol = N)
  for (t in 2:(n + burn)) {
    xi[t, ] = xi[t-1, ] %*% t(A) + e[t, ]
  }
  xi <- xi[1:n + burn, ]
  
  s <- sapply(1:N, function(i){sd(xi[, i]) / sd(chi[, i])})
  y <- chi * matrix(s, nrow = n, ncol = N, byrow = TRUE) + xi
  out <- list("x" = y, "A1" = A, "lambda_max" = max_ev_A, "D" = D, "Lambda" = Lambda, "u"=e)
}

sim_DGP <- function(n, N, type = 1, nburn = 50, mu = 0, prop = 1, ...) {
  if (type == 0) {
    sim <- sim_DGP_FNETS(n = n, N = N, burn = nburn, ...)
  } else if (type > 0) {
    sim <- sim_DGP_VAR(n = n, N = N, type = type, nburn = nburn, ...)
  }
  means <- rep(0, N)
  means[1:ceiling(prop * N)] <- mu
  sim$x <- sim$x + matrix(means, ncol = N, nrow = n, byrow = TRUE)
  return(list("x"=sim$x, "A"=sim$A, "u"=sim$u))
}
  