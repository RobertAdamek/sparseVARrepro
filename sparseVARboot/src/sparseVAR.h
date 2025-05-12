#include <RcppArmadillo.h>

#ifndef sparseVAR_H
#define sparseVAR_H

struct VAR_out {
  arma::mat coef;
  arma::mat resid;
  /////////////remove
  double lambda;
  arma::vec lambdas;
  ///////////////////
};

struct VAR_select_out {
  arma::mat coef;
  arma::mat resid;
  double lambda;
  arma::vec lambdas;
};

VAR_out sparseVAR(arma::mat Y, const int& p, const bool& trim, const int& pen, const double& nbr_lambdas,
                  const double& lambda_ratio, const double& eps, const int& selection,
                  const bool& pen_own, const bool& only_lag1, const double& c, 
                  const unsigned int K = 15, double improvement_thresh = 0.01, unsigned int Nsim=1000, const double& alpha = 0.05);

arma::mat lag_matrix(const arma::mat& x, const int& p, const bool& trim);

double mySquare(double x);

#endif

