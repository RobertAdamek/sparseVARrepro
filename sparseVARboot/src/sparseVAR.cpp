#include <RcppArmadillo.h>
#include "sparseVAR.h"

using namespace arma;
using namespace std;

// Functions for L1 penalization //


double ST1a(const double& z, const double& gam){
  
  if(z>0 && gam<fabs(z)) return(z-gam);
  if(z<0 && gam<fabs(z)) return(z+gam);
  if(gam>=fabs(z)) return(0);
  else return(0);
}

arma::colvec ST3a(const arma::colvec& z, const double& gam, const int& index_own, const int& p, const bool& pen_own,
                  const bool& only_lag1){
  
  // pen_own: boolean : TRUE if penalize own lags, FALSE to leave own lags unpenalized
  // only_lag1 : boolean, only relevant if pen_own = false : TRUE if only first lag should be unpenalized, FALSE all own lags should be unpenalized
  
  int n=z.size(); // number of time series times number of lags
  int k = n/p; // number of time series
  arma::colvec var_index = zeros(k); 
  var_index(index_own) = 1; //indicate 1 at spot of own time series
  
  arma::colvec index_own_lags = var_index; 
  if(only_lag1){
    index_own_lags = join_vert(index_own_lags, zeros((p-1)*k)); 
  }else{
    
    for(int j=0; j<p-1; ++j){
      index_own_lags = join_vert(index_own_lags, var_index);
    }
  }
  
  
  
  
  arma::colvec z1(n);
  for( int i=0; i<n;++i){
    double z11=z(i);
    
    if((!pen_own) & (index_own_lags(i)==TRUE)){
      z1(i) = ST1a(z11, 0);
    }else{
      z1(i) = ST1a(z11, gam);
    }
  }
  return(z1);
}

arma::mat FistaLV(const arma::mat& Y, const arma::mat& Z, arma::mat& B, const double& gam, 
                  const double& eps, const double& tk, const int& k, const int& p, const bool& pen_own, const bool& only_lag1){
  B=trans(B);
  arma::colvec B1=B.col(0);
  
  double j = 1;
  
  for( int i =0; i<k; ++i){
    B1=B.col(i);
    arma::colvec BOLD=B.col(i);
    arma::colvec BOLDOLD=BOLD;
    double thresh=10*eps;
    j=1;
    
    while(thresh>eps){
      
      arma::colvec v=BOLD+((j-2)/(j+1))*(BOLD-BOLDOLD);
      
      B1=ST3a(vectorise(v)+tk*vectorise((trans(Y.col(i))-trans(v)*Z)*trans(Z)),gam*tk, i, p, pen_own, only_lag1);
      thresh=max(abs(B1-v));
      BOLDOLD=BOLD;
      BOLD=B1;
      j+=1;
      if(j>10000){
        break;
      }
    }
    
    B.col(i)=B1;
  }
  
  B=trans(B);
  return(B);
}

arma::cube gamloopFista2(const arma::cube& bcube, const arma::mat& Y,const arma::mat& Z, 
                         const arma::colvec& gammgrid, const double& eps,
                         const arma::colvec& YMean2, const arma::colvec& ZMean2, 
                         arma::mat& B1, const int& k, const int& p, const double& tk,
                         const bool& pen_own, const bool& only_lag1) {
  
  arma::mat b2=B1;
  arma::mat B1F2=B1;
  int ngridpts = bcube.n_slices;
  arma::cube bcube2(k, k*p + 1, ngridpts);
  bcube2.fill(0);
  arma::colvec nu=zeros<colvec>(k);
  
  double gam =0;
  
  int i;
  //loop through candidate lambda values
  for (i=0; i<ngridpts;++i) {
    gam=gammgrid[i];
    arma::mat B1F2 = bcube.slice(i);
    B1 = FistaLV(Y, Z, B1F2, gam, eps, tk, k, p, pen_own, only_lag1);
    nu = YMean2 - B1 *ZMean2;
    bcube2.slice(i) = mat(join_horiz(nu, B1));
  }
  
  return(bcube2);
}

arma::cube lassoVARFistcpp(const arma::cube& beta, const arma::mat& trainY, const arma::mat& trainZ, const arma::colvec& lambda,
                           const double& tol, const int& p, const bool& pen_own, const bool& only_lag1, const double& tk){
  
  int n = trainY.n_rows;
  int k = trainY.n_cols;
  int g = beta.n_slices;
  arma::mat YMean = mean(trainY);
  arma::mat ZMean = mean(trainZ.t());
  
  arma::mat Y = zeros(n, k);
  arma::mat trainZt = trainZ.t();
  arma::mat Z = zeros(n, k*p);
  for (int i = 0; i < n; i++) {
    Y.row(i) = trainY.row(i)  - YMean ;
    Z.row(i) = trainZt.row(i)  - ZMean ;
  }
  Z = Z.t();
  
  // // Step size
  // vec eigval;
  // arma::mat eigvec;
  // const arma::mat Zt=Z*trans(Z);
  // eig_sym(eigval, eigvec, Zt);
  // double tk = 1/max(eigval);
  
  // Rcpp::Rcout << "tk computed " << tk << std::endl;
  
  arma::mat betaini1 = beta.subcube(0, 1, 0, k-1, k*p + 1 - 1, 0);
  arma::cube betaini = beta.subcube(0, 1, 0, k-1, k*p + 1 - 1, g-1);
  
  arma::cube betaout = gamloopFista2(betaini, Y, Z, lambda, tol, YMean.t(), ZMean.t(), betaini1, k, p, tk, pen_own, only_lag1);
  
  return(betaout);
}


// Functions for HLag penalization //
uvec bbsubs(const int& j, const int& k, const int& p)
{
  uvec bb(p);
  bb(0)=j;
  for(int i=1;i<p;++i){
    bb(i)=j+k*(i);
  }
  
  return(bb);
}

uvec vsubscppelem(const int& p, const int& pmax){
  uvec vs(pmax-p+1);
  for(int i=pmax;i>=p;--i){
    vs(i-p)=i-1;
  }
  return(vs);
}

rowvec proxcppelem(const arma::colvec& v2, const int& L, const double& lambda,
                   const uvec& res1, const arma::colvec& w){
  arma::colvec r =v2;
  for(int i=(L-1); i>=0;--i){
    
    uvec res=vsubscppelem(i+1,L);
    
    if(norm(r(res)/(lambda*w(i)),"fro")<1+1e-8){
      r(res)=zeros(res.n_elem);
    }else{
      r(res)=r(res)-lambda*w(i)*r(res)/(norm(r(res),"fro"));
    }
  }
  return(trans(r));
}

arma::rowvec prox2(const arma::colvec& v, const double& lambda, const int& k, const int& p,
                   const uvec& res1, const arma::colvec& w)
{
  arma::rowvec v2(v.n_elem);
  arma::rowvec v3(p);
  for(int i=0;i<k;++i)
  {
    uvec bb=bbsubs(i,k,p);
    arma::colvec v1=v(bb);
    v3=proxcppelem(v1,p,lambda,res1,w);
    v2(bb)=v3;
  }
  return(v2);
}


// double norm2(Rcpp::NumericVector x){
//   arma::vec xx = x;
//   double g=arma::norm(xx,2);
//   return (Rcpp::as<double>(Rcpp::wrap(g)));
// }

uvec ind(int n2,int m){
  
  std::vector<int> subs;
  
  for(int i =0 ; i<n2;++i){
    
    subs.push_back(i);
    
  }
  
  subs.erase(subs.begin()+m);
  
  return(conv_to<uvec>::from(subs));
}


arma::mat FistaElem(const arma::mat& Y,const arma::mat& Z, const arma::mat& phi, const int p,const int k,double lambda, const double eps,const double tk){
  double j=1;
  arma::mat phiFin=phi;
  arma::rowvec phiR=phi.row(0);
  arma::rowvec phiOLD=phiR;
  arma::rowvec phiOLDOLD=phiOLD;
  arma::rowvec v=phiOLD;
  uvec res1=ind(p,0);
  arma::colvec w(p);
  w.ones();
  
  for(int i=0;i<k;++i)
  {
    j=1;
    double thresh=10*eps;
    phiR=phi.row(i);
    phiOLD=phiR;
    phiOLDOLD=phiOLD;
    v=phiR;
    while(thresh>eps)
    {
      v=phiOLD+((j-2)/(j+1))*(phiOLD-phiOLDOLD);
      
      phiR=prox2(vectorise(v)+tk*vectorise((trans(Y.col(i))-v*Z)*trans(Z)),tk*lambda,k,p,res1,w);
      
      thresh=max(abs(phiR-v));
      phiOLDOLD=phiOLD;
      phiOLD=phiR;
      j+=1;
      
    }
    phiFin.row(i)=phiR;
  }
  return(phiFin);
}

arma::cube gamloopElem2(const arma::cube &bcube, const arma::mat& Y,const arma::mat& Z, const arma::colvec& gammgrid, const double& eps,
                        const arma::colvec& YMean2, const arma::colvec& ZMean2, arma::mat &B1, const int& k, const int& p, const double& tk,
                        const int& flag_restart_opt = 1){
  
  arma::mat B1F2 = B1;
  int ngridpts = bcube.n_slices;
  arma::cube bcube2(k, k*p + 1, ngridpts);
  bcube2.fill(0);
  colvec nu=zeros<colvec>(k);
  
  int i;
  B1 = bcube.slice(0);
  for (i=0; i<ngridpts;++i) {
    B1F2 = (flag_restart_opt == 1) ? B1 : bcube.slice(i);
    B1 = FistaElem(Y, Z, B1F2, p, k, gammgrid[i], eps, tk);
    nu = YMean2 - B1*ZMean2;
    bcube2.slice(i) = mat(join_horiz(nu, B1));
  }
  
  return(bcube2);
}

arma::cube HVARElemAlgcpp(const arma::cube &beta, const arma::mat& trainY, const arma::mat& trainZ, const arma::colvec& lambda,
                          const double& tol, const int& p, const double& tk, const int& flag_restart_opt = 0){
  // Prelimaries
  int n = trainY.n_rows;
  int k = trainY.n_cols;
  int g = beta.n_slices;
  arma::mat YMean = mean(trainY);
  arma::mat ZMean = mean(trainZ.t());
  
  arma::mat Y = zeros(n, k);
  arma::mat trainZt = trainZ.t();
  arma::mat Z = zeros(n, k*p);
  for (int i = 0; i < n; i++) {
    Y.row(i) = trainY.row(i)  - YMean ;
    Z.row(i) = trainZt.row(i)  - ZMean ;
  }
  Z = Z.t();
  
  // // Step size
  // vec eigval;
  // arma::mat eigvec;
  // const arma::mat Zt=Z*trans(Z);
  // eig_sym(eigval, eigvec, Zt);
  // double tk = 1/max(eigval);
  // 
  // Rcpp::Rcout << "tk computed " << tk << std::endl;
  
  arma::mat betaini1 = beta.subcube(0, 1, 0, k-1, k*p + 1 - 1, 0);
  const arma::cube betaini = beta.subcube(0, 1, 0, k-1, k*p + 1 - 1, g-1);
  arma::cube betaout = gamloopElem2(betaini, Y, Z, lambda, tol, YMean.t(), ZMean.t(), betaini1, k, p, tk, flag_restart_opt);
  
  return(betaout);
}


bool moveup_LGSearch_cpp(const arma::mat &param){
  int n = param.n_rows;
  int k = param.n_cols;
  
  for (int i=0; i<n; i++){
    for (int j=0; j<k; j++){
      if (param(i, j) != 0) return(false);
    }
  }
  return(true);
}


double LGSearch_cpp(const double& gstart, const arma::mat &Y, const arma::mat &Z, 
                    arma::cube beta, const int& estim, const int& k, const int& p, const bool& pen_own, 
                    const bool& only_lag1, const double& tk) {
  double lambdah = gstart;
  double lambdal = 0.0;
  arma::mat param;
  
  while (std::abs(lambdah - lambdal) > 0.00001){
    double lambda = (lambdah + lambdal)/2;
    arma::colvec lvec(1);
    lvec.fill(lambda);
    if (estim == 1){
      beta = lassoVARFistcpp(beta, Y, Z, lvec, 0.0001, p, pen_own, only_lag1, tk);
      param = beta.slice(0).cols(1, k*p);
    }
    else if (estim == 2){
      beta = HVARElemAlgcpp(beta, Y, Z, lvec, 0.0001, p, tk);
      param = beta.slice(0).cols(1, k*p);
    }
    
    bool move_up = moveup_LGSearch_cpp(param);
    if (move_up) lambdah = lambda;
    else lambdal = lambda;
    
  }
  
  return(lambdah);
}

double VAR_ic(const arma::mat& res, const arma::mat& coef, const int& selection){
  // selection : 1 bic; 2 aic; 3 hq
  
  arma::vec coef_nz = nonzeros(coef);
  double df = coef_nz.size();
  double tt = res.n_rows;
  arma::mat Omega = cov(res);
  double ic;
  
  if(selection==1){
    ic = log(det(Omega)) + (log(tt)/tt) * df;
  }
  if(selection==2){
    ic = log(det(Omega)) + (2/tt) * df;
  }
  if(selection==3){
    ic = log(det(Omega)) + (2*log(log(tt))/tt) * df;
  }
  return(ic);
}

VAR_out HVAR(const arma::mat& fullY, const arma::mat& fullZ, const int& p, const int& k, 
             const arma::colvec& lambda, const double& eps, const int& pen, const bool& pen_own, 
             const bool& only_lag1, const double& tk){
  // fullY matrix with tt rows and k columns
  // fullZ matrix with kp rows and tt columns
  // pen 1: L1; 2: HLag
  
  arma::cube betas(k, k*p+1, 1, fill::zeros);
  arma::cube Phis(k, k*p+1, 1);
  
  if(pen==1){
    Phis = lassoVARFistcpp(betas, fullY, fullZ, lambda, eps, p, pen_own, only_lag1, tk);
  }
  
  if(pen==2){
    Phis = HVARElemAlgcpp(betas, fullY, fullZ, lambda, eps, p, tk);
  }
  
  const arma::mat fullZi = join_horiz(ones(fullZ.n_cols), fullZ.t());
  VAR_out out;
  out.resid = (fullY.t() - Phis.slice(0)*fullZi.t()).t();
  arma::mat coef = Phis.slice(0);
  out.coef = coef.submat(0, 1, k-1,p*k);
  
  return out;
  // return Rcpp::List::create(
  //   Rcpp::Named("resid") = resid,
  //   Rcpp::Named("coef") = coef.submat(0, 1, k-1,p*k));
  
  // output slot coef has (k) rows and (kp) columns!
}


arma::vec seq_default(long double from, long double to, long unsigned int length_out){
  return arma::linspace(from, to, length_out);
}

arma::vec LambdaGridE(double gran1, double gran2, arma::mat Y, arma::mat Z, int pen, int p, int k,
                      const bool& pen_own, const bool& only_lag1, const double& tk){
  // gran1 : ratio lambda_max/lambda_min
  // gran2 : number of lambdas
  // pen 1: L1; 2: HLag
  
  arma::vec lambda_starts(k);
  arma::mat crossprod;
  double lambda_start;
  
  if(pen==1){
    crossprod = Y.t()*Z.t();
    lambda_start = crossprod.max();
  }
  
  if(pen==2){
    for(int i=0; i < k;++i){
      lambda_starts(i) = arma::norm(Z *Y.col(i), 2);
    }
    lambda_start = lambda_starts.max();
  }
  
  arma::cube betas(k, k*p+1, 1, fill::zeros);
  
  lambda_start = LGSearch_cpp(lambda_start, Y, Z, betas, pen, k, p, pen_own, only_lag1, tk);
  
  arma::vec lambda_grid = seq_default(log(lambda_start), log(lambda_start/gran1), gran2);
  lambda_grid = exp(lambda_grid);
  
  return(lambda_grid);
}


arma::mat lag_matrix(const arma::mat& x, const int& p, const bool& trim) {//step 3 in the bootstrap algorithm: used with trim=false to make the regressor matrix
  const int n = x.n_rows;
  const int k = x.n_cols;
  arma::mat lag_x = arma::mat(n, k * p, fill::zeros); 
  for (int j = 0; j < p; j++) {
    lag_x.submat(j + 1, k * j, n - 1, k * (j + 1) - 1) = x.rows(0, n - j - 2);
  }
  return lag_x.rows(p * trim, n - 1);
}

double Andrews91_truncation(const arma::mat& What, const unsigned int& T_, const unsigned int& h){
  arma::vec rhos(h), variances(h);
  arma::vec y_T(T_-1), y_Tm1(T_-1), constant(T_-1,fill::ones), residual(T_-1);
  arma::mat X_T(T_-1,2);
  arma::vec beta(2);
  for(unsigned int i=0; i<h; i++){
    y_T=What.submat(1,i,T_-1,i);
    y_Tm1=What.submat(0,i,T_-2,i);
    X_T=join_horiz(constant,y_Tm1);
    beta=inv(X_T.t()*X_T)*X_T.t()*y_T;
    residual=y_T-X_T*beta;
    rhos(i)=beta(1);
    variances(i)=as_scalar(residual.t()*residual)/double(double(T_)-2.0);
  }
  double numerator=0, denominator=0;
  for(unsigned int i=0; i<h; i++){
    numerator+=4*pow(rhos(i),2)*pow(variances(i),2)/double(pow(1-rhos(i),6)*pow(1+rhos(i),2));
    denominator+=pow(variances(i),2)/double(pow(1-rhos(i),4));
  }
  double alphahat1=numerator/double(denominator);
  double S_T=1.1447*pow(alphahat1*double(T_),double(1.0/double(3.0)));
  return S_T;
}

arma::mat LRVestimator(const arma::mat& U, const arma::mat& X, const double& LRVtrunc, const double& T_multiplier){
  // U : T_ by d matrix of VAR residuals for VAR with d time series
  // X : T_ by N matrix of VAR predictors, for a VAR with d time series and p lags, N=dp
  
  const unsigned int& d = U.n_cols; // number of time series
  const unsigned int& N = X.n_cols; // N is number of predictor variables (pxd in VAR)
  const unsigned int& T_ = X.n_rows; // number of time points
  const unsigned int& Nd = N*d; // dimension of the long-run variance
  
  arma::mat Omegahat(Nd, Nd,fill::zeros);
  arma::mat What(T_, Nd);
  arma::mat Xi_ell(Nd, Nd);
  unsigned int counter = 0;
  for(unsigned int ui=0; ui < d; ui++){
    for(unsigned int xj=0; xj < N; xj++){
      What.col(counter)= U.col(ui)%X.col(xj);
      counter = counter + 1;
    }
  }
  double Q_T_double;
  int Q_T;
  if(LRVtrunc==0 && T_multiplier==0){//If both the T_multiplier and LRVtrunc are 0, do a data-driven choice of tuning parameter
    Q_T_double=std::ceil(Andrews91_truncation(What, T_, Nd));
    Q_T= (int) Q_T_double;
  }else{
    Q_T_double=std::ceil(pow(T_multiplier*double(T_),LRVtrunc));
    Q_T= (int) Q_T_double;
  }
  if(Q_T_double>double(T_)/2.0){
    Q_T_double=std::ceil(double(T_)/2.0);
    Q_T= (int) Q_T_double;
  }
  if(Q_T<=1){ //If Q_T is zero or negative, Omegahat is Xi(0)
    Omegahat=(1.0/double(T_))*What.t()*What;
  }else{ //Otherwise
    Omegahat=(1.0/double(T_))*What.t()*What;
    for(int ell=1; ell<=Q_T-1; ell++){
      Xi_ell=(1.0/double(T_-ell))*(What.rows(ell,T_-1)).t()*What.rows(0,T_-ell-1);
      Omegahat+=(1-double(ell)/double(Q_T))*(Xi_ell+Xi_ell.t());
    }
  }
  return(Omegahat);
}

VAR_select_out selectPI(const arma::vec lambda_grid, const arma::mat Y, const arma::mat VAR_lags, const int& p, const int& k, 
                        const double& eps, const int& pen, const bool& pen_own, const bool& only_lag1, const double& c, 
                        const unsigned int K, double improvement_thresh, unsigned int B, const double& alpha, const double& tk){
  // c constant for plug-in approach. Consider 0.4/0.8
  // K is the number of iterations of the plug-in approach
  // improvement_thresh if the % change in lambda is less than this, stop iterating
  // B number of simulations used to estimate the quantiles of the Gaussian maximum
  // alpha: alpha quantile of the Gaussian maximum to be considered
  
  // Plug-in approach to select one lambda for the VAR
  
  arma::mat LRCovariance; 
  double lambda_old = 0; // need double to check convergence
  lambda_old = lambda_grid(0);
  arma::colvec lambda_old_vec; // need vec as input to HVAR
  lambda_old_vec = lambda_grid(0);
  VAR_out HVAR_constant;
  HVAR_constant =  HVAR(Y, VAR_lags.t(), p, k, lambda_old_vec, eps, pen, true, false, tk); // get residuals from constant only model
  arma::mat uhat = HVAR_constant.resid; 
  
  const unsigned int& d = uhat.n_cols; // number of time series
  const unsigned int& N = VAR_lags.n_cols; // N is number of predictor variables (pxd in VAR)
  const unsigned int& T_ = VAR_lags.n_rows; // number of time points
  const unsigned int& Nd = N*d; // dimension of the long-run variance
  
  double lambda;
  arma::colvec lambda_vec(1);
  lambda_vec(0) = 0;
  VAR_out HVAR_fit;
  arma::vec Gmax(B);
  arma::mat Gs;
  arma::vec Gmeans(Nd); 
  double cutoff;
  arma::vec eigval(Nd);
  arma::mat eigvec(Nd,Nd);
  arma::mat sqrt_cov(Nd,Nd);
  arma::mat sqrt_diag_eigval(Nd,Nd,fill::zeros);
  arma::mat random_gaussians=randn(Nd,B);
  
  //unsigned int iteration=K;
  
  //loop where we iterate to find the best lambda
  for(unsigned int ik=0; ik< K; ik++){ 
    
    // Rcpp::Rcout << "iteration " << ik << std::endl;
    
    LRCovariance = LRVestimator(uhat, VAR_lags, 0, 0); // returns an Nd by Nd matrix

    eig_sym(eigval, eigvec, LRCovariance);

    for(unsigned int i=0; i<Nd; i++){
      sqrt_diag_eigval(i,i)=sqrt(max(0.0,eigval(i))); //negative eigenvalues replaced by 0
    }

    sqrt_cov=eigvec*sqrt_diag_eigval;
    
    Gs=sqrt_cov*random_gaussians; //generate the correlated gaussians

    for(unsigned int b=0; b<B; b++){
      Gmax(b)=max(abs(Gs.col(b)));
    }

    Gmax=sort(Gmax);
    cutoff=Gmax(B*(1-alpha)); //1-alpha quantile of the randomly generated data
    lambda=c*cutoff*double(sqrt(double(T_))); ///not yet multiplying by 4 here
    lambda_vec(0) = lambda;
    
    // Rcpp::Rcout << "lambda " << lambda << std::endl;
    
    if(abs(lambda-lambda_old)/lambda_old<improvement_thresh){ //Check if the improvement is big enough
      //iteration=k;
      ik=K; //no more loops after this
    }
    
    HVAR_fit =  HVAR(Y, VAR_lags.t(), p, k, lambda_vec, eps, pen, pen_own, only_lag1, tk); // get residuals from sparse VAR
    uhat = HVAR_fit.resid;
    lambda_old = lambda;
    lambda_old_vec = lambda_vec;
    
  }
  
  const arma::mat coef = HVAR_fit.coef;
  const arma::mat resid = HVAR_fit.resid;
  
  VAR_select_out out;
  out.coef = coef.t();
  out.resid = resid;
  out.lambda = lambda;
  out.lambdas = lambda_grid;
  
  return out; 
}

// added by Robert
double stdnormal_cdf(double u)
{
  const double a[5] = {
    1.161110663653770e-002,3.951404679838207e-001,2.846603853776254e+001,
    1.887426188426510e+002,3.209377589138469e+003
  };
  const double b[5] = {
    1.767766952966369e-001,8.344316438579620e+000,1.725514762600375e+002,
    1.813893686502485e+003,8.044716608901563e+003
  };
  const double c[9] = {
    2.15311535474403846e-8,5.64188496988670089e-1,8.88314979438837594e00,
    6.61191906371416295e01,2.98635138197400131e02,8.81952221241769090e02,
    1.71204761263407058e03,2.05107837782607147e03,1.23033935479799725E03
  };
  const double d[9] = {
    1.00000000000000000e00,1.57449261107098347e01,1.17693950891312499e02,
    5.37181101862009858e02,1.62138957456669019e03,3.29079923573345963e03,
    4.36261909014324716e03,3.43936767414372164e03,1.23033935480374942e03
  };
  const double p[6] = {
    1.63153871373020978e-2,3.05326634961232344e-1,3.60344899949804439e-1,
    1.25781726111229246e-1,1.60837851487422766e-2,6.58749161529837803e-4
  };
  const double q[6] = {
    1.00000000000000000e00,2.56852019228982242e00,1.87295284992346047e00,
    5.27905102951428412e-1,6.05183413124413191e-2,2.33520497626869185e-3
  };
  double y, z;
  
  if (std::isnan(u))
    return std::nan("");
  if (!std::isfinite(u))
    return (u < 0 ? 0.0 : 1.0);
  y = fabs(u);
  if (y <= 0.46875*M_SQRT2) {
    /* evaluate erf() for |u| <= sqrt(2)*0.46875 */
    z = y*y;
    y = u*((((a[0]*z+a[1])*z+a[2])*z+a[3])*z+a[4])
      /((((b[0]*z+b[1])*z+b[2])*z+b[3])*z+b[4]);
    return 0.5+y;
  }
  z = exp(-y*y/2)/2;
  if (y <= 4.0) {
    /* evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0 */
    y = y/M_SQRT2;
    y =
      ((((((((c[0]*y+c[1])*y+c[2])*y+c[3])*y+c[4])*y+c[5])*y+c[6])*y+c[7])*y+c[8])
      
      
      /((((((((d[0]*y+d[1])*y+d[2])*y+d[3])*y+d[4])*y+d[5])*y+d[6])*y+d[7])*y+d[8]);
    
    y = z*y;
  } else {
    /* evaluate erfc() for |u| > sqrt(2)*4.0 */
    z = z*M_SQRT2/y;
    y = 2/(y*y);
    y = y*(((((p[0]*y+p[1])*y+p[2])*y+p[3])*y+p[4])*y+p[5])
      /(((((q[0]*y+q[1])*y+q[2])*y+q[3])*y+q[4])*y+q[5]);
    y = z*(1.0/sqrt(M_PI)-y);
  }
  return (u < 0.0 ? y : 1-y);
};

double stdnormal_inv(double p)
{
  const double a[6] = {
    -3.969683028665376e+01,  2.209460984245205e+02,
    -2.759285104469687e+02,  1.383577518672690e+02,
    -3.066479806614716e+01,  2.506628277459239e+00
  };
  const double b[5] = {
    -5.447609879822406e+01,  1.615858368580409e+02,
    -1.556989798598866e+02,  6.680131188771972e+01,
    -1.328068155288572e+01
  };
  const double c[6] = {
    -7.784894002430293e-03, -3.223964580411365e-01,
    -2.400758277161838e+00, -2.549732539343734e+00,
    4.374664141464968e+00,  2.938163982698783e+00
  };
  const double d[4] = {
    7.784695709041462e-03,  3.224671290700398e-01,
    2.445134137142996e+00,  3.754408661907416e+00
  };
  
  double q, t, u;
  
  if (std::isnan(p) || p > 1.0 || p < 0.0)
    return std::nan("");
  if (p == 0.0)
    return std::numeric_limits<double>::infinity();
  if (p == 1.0)
    return std::numeric_limits<double>::infinity();
  q = std::min(p,1-p);
  if (q > 0.02425) {
    /* Rational approximation for central region. */
    u = q-0.5;
    t = u*u;
    u = u*(((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4])*t+a[5])
      /(((((b[0]*t+b[1])*t+b[2])*t+b[3])*t+b[4])*t+1);
  } else {
    /* Rational approximation for tail region. */
    t = sqrt(-2*log(q));
    u = (((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5])
      /((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1);
  }
  /* The relative error of the approximation has absolute value less
   than 1.15e-9.  One iteration of Halley's rational method (third
   order) gives full machine precision... */
  t = stdnormal_cdf(u)-q;    /* error */
  t = t*M_SQRT2*sqrt(M_PI)*exp(u*u/2);   /* f(u)/df(u) */
  u = u-t/(1+u*t/2);     /* Halley's method */
  
  return (p > 0.5 ? -u : u);
};

struct lasso_output{
  unsigned int N, T_, gridsize;
  arma::vec grid, y;
  arma::mat betahats, X;
  int opt_type; //1="naive", 2="covariance", 3="adaptive"
};

struct partial_lasso_output{
  bool partial;
  unsigned int N, T_, gridsize, h;
  arma::vec grid, y;
  arma::uvec H, minusH;
  arma::mat betahats, betahats_1, betahats_2, X, X_1, X_2;
  int opt_type; //1="naive", 2="covariance", 3="adaptive"
};

double soft_threshold(const double& z, const double& gamma){
  double ret;
  if(std::abs(z)<=gamma){
    ret=0;
  }else if(z-gamma>0){
    ret=z-gamma;
  }else{
    ret=z+gamma;
  }
  return ret;
}

arma::mat coordinate_descent_naive(const arma::mat& X, const arma::colvec& y, const arma::vec& grid, const double& opt_threshold,
                                   const unsigned int& N, const unsigned int& T_, const unsigned int& gridsize){
  arma::mat betahats(N,gridsize);
  arma::vec betahat(N,fill::zeros);
  arma::vec betahat_old(N,fill::zeros);
  unsigned int g,j,k;
  arma::uvec Active_set; Active_set.reset(); //see if this is necessary
  arma::uvec temp_Active_set;
  unsigned int active_length=0, nonactive_length=N;
  arma::uvec Nonactive_set=linspace<arma::uvec>(0,N-1,N);
  double change=opt_threshold+1;
  arma::vec x_crossprods(N,fill::zeros);
  arma::vec x_j(T_);
  for(j=0;j<N;j++){
    x_j=X.col(j);
    x_crossprods(j)=as_scalar(x_j.t()*x_j);
  }
  arma::vec full_res=y;
  arma::vec partial_res=full_res;
  //loop through lambda grid
  for(g=0;g<gridsize;g++){
    //run through non-active set
    //Nonactive_set=linspace<arma::uvec>(0,N-1,N); Nonactive_set.shed_rows(Active_set);
    nonactive_length=Nonactive_set.n_elem;
    for(k=0; k<nonactive_length;k++){
      j=Nonactive_set(k);
      x_j=X.col(j);
      partial_res=full_res;//+X.col(j)*betahat(j);
      betahat(j)=soft_threshold(as_scalar(x_j.t()*partial_res)/x_crossprods(j), T_*grid(g)/x_crossprods(j));
      if(betahat(j)!=0){ //new variable becomes nonzero
        active_length++;
        temp_Active_set=arma::uvec(active_length);
        temp_Active_set.head_rows(active_length-1)=Active_set;
        temp_Active_set(active_length-1)=j;
        Active_set=sort(temp_Active_set); //sorting seems to make it faster for some reason
        full_res=partial_res-X.col(j)*betahat(j);// only need to update if this is nonzero
      }
    }
    Nonactive_set=linspace<arma::uvec>(0,N-1,N); Nonactive_set.shed_rows(Active_set); //update nonactive set
    change=opt_threshold+1;
    while(change>opt_threshold){
      //for loop through parameters of the active set
      for(k=0;k<active_length;k++){
        j=Active_set(k);
        x_j=X.col(j);
        partial_res=full_res+X.col(j)*betahat(j);
        betahat(j)=soft_threshold(as_scalar(x_j.t()*partial_res)/x_crossprods(j), T_*grid(g)/x_crossprods(j));
        full_res=partial_res-X.col(j)*betahat(j);
      }
      //change=sqrt(as_scalar((betahat-betahat_old).t()*(betahat-betahat_old)))/double(N);
      change=max(abs(betahat-betahat_old));
      betahat_old=betahat;
    }
    betahats.col(g)=betahat;
  }
  return betahats;
}

arma::mat coordinate_descent_covariance(const arma::mat& X, const arma::colvec& y, const arma::vec& grid, const double& opt_threshold,
                                        const unsigned int& N, const unsigned int& T_, const unsigned int& gridsize){
  arma::mat betahat_mat(N,gridsize);
  arma::vec betahat(N,fill::zeros);
  arma::vec betahat_old(N,fill::zeros);
  unsigned int g,j,k,l,a,b;
  arma::uvec Active_set;
  arma::uvec temp_Active_set;
  unsigned int active_length=0, nonactive_length=N;
  arma::uvec Nonactive_set=linspace<arma::uvec>(0,N-1,N);
  double change=opt_threshold+1;
  arma::mat x_crossprods(N,N,fill::zeros);
  arma::vec x_j(T_);
  for(j=0;j<N;j++){
    x_j=X.col(j);
    x_crossprods(j,j)=as_scalar(x_j.t()*x_j);
  }
  arma::vec xy(N,1);
  for(j=0;j<N;j++){
    x_j=X.col(j);
    xy(j)=as_scalar(x_j.t()*y);
  }
  //loop through lambda grid
  for(g=0;g<gridsize;g++){
    //one run through non-active set
    nonactive_length=Nonactive_set.n_elem;
    for(k=0; k<nonactive_length;k++){
      j=Nonactive_set(k);
      double x_jx_kbeta_k=0;
      for(l=0;l<N;l++){
        if(l!=j && betahat(l)!=0){
          x_jx_kbeta_k+=x_crossprods(j,l)*betahat(l);
        }
      }
      betahat(j)=soft_threshold((xy(j)-x_jx_kbeta_k)/x_crossprods(j,j), T_*grid(g)/x_crossprods(j,j));
      if(betahat(j)!=0){ //new variable becomes nonzero
        active_length++;
        x_j=X.col(j);
        for(l=0; l<N;l++){ //calculate crossproducts
          x_crossprods(j,l)=x_crossprods(l,j)=as_scalar(x_j.t()*X.col(l));
        }
        temp_Active_set=arma::uvec(active_length);
        temp_Active_set.head_rows(active_length-1)=Active_set;
        temp_Active_set(active_length-1)=j;
        Active_set=sort(temp_Active_set); //sorting seems to make it faster for some reason
      }
      
    }
    Nonactive_set=linspace<arma::uvec>(0,N-1,N); Nonactive_set.shed_rows(Active_set); //update nonactive set
    change=opt_threshold+1;
    //optimize
    while(change>opt_threshold){
      //loop through parameters of the active set
      for(k=0;k<active_length;k++){
        j=Active_set(k);
        x_j=X.col(j);
        double x_jx_kbeta_k=0;
        for(a=0; a<active_length;a++){
          b=Active_set(a);
          if(b!=j){
            x_jx_kbeta_k+=x_crossprods(j,b)*betahat(b);
          }
        }
        betahat(j)=soft_threshold((xy(j)-x_jx_kbeta_k)/x_crossprods(j,j), T_*grid(g)/x_crossprods(j,j));
      }
      //change=sqrt(as_scalar((betahat-betahat_old).t()*(betahat-betahat_old)))/N;
      change=max(abs(betahat-betahat_old));
      betahat_old=betahat;
    }
    betahat_mat.col(g)=betahat;
  }
  return betahat_mat;
}

lasso_output lasso(const arma::mat& X, const arma::colvec& y, const arma::vec& grid,
                   const double& opt_threshold, const int& opt_type){
  unsigned int T_=X.n_rows;
  unsigned int N=X.n_cols;
  unsigned int gridsize=grid.n_elem;
  arma::mat betahats(N,gridsize);
  switch(opt_type) {
  case 1: //"naive"
    betahats=coordinate_descent_naive(X, y, grid, opt_threshold,
                                      N, T_, gridsize);
    break;
  case 2: //"covariance"
    betahats=coordinate_descent_covariance(X, y, grid, opt_threshold,
                                           N, T_, gridsize);
    break;
  case 3: //"adaptive"
    if(N>T_){
      betahats=coordinate_descent_naive(X, y, grid, opt_threshold,
                                        N, T_, gridsize);
    }
    else{
      betahats=coordinate_descent_covariance(X, y, grid, opt_threshold,
                                             N, T_, gridsize);
    }
    break;
  default:
    //warning("Warning: Invalid opt_type, choosing type 3");
    if(N>T_){
      betahats=coordinate_descent_naive(X, y, grid, opt_threshold,
                                        N, T_, gridsize);
    }
    else{
      betahats=coordinate_descent_covariance(X, y, grid, opt_threshold,
                                             N, T_, gridsize);
    }
  }
  lasso_output ret;
  ret.N=N;
  ret.T_=T_;
  ret.gridsize=gridsize;
  ret.grid=grid;
  ret.y=y;
  ret.betahats=betahats;
  ret.X=X;
  ret.opt_type=opt_type;
  return(ret);
}

lasso_output lasso_weighted(const arma::mat& X, const arma::colvec& y, const arma::vec& grid, const arma::vec& weights, 
                            const double& opt_threshold, const int& opt_type){
  // this is based on https://stats.stackexchange.com/questions/397986/how-to-solve-an-adaptive-lasso-model
  unsigned int N=X.n_cols;
  arma::mat X_weighted=X; 
  for(unsigned int i=0; i<N; i++){
    X_weighted.col(i)=X.col(i)*(1/weights(i));
  }
  lasso_output L_weighted=lasso(X_weighted, y, grid, opt_threshold, opt_type);
  for(unsigned int i=0; i<N; i++){
    L_weighted.betahats.row(i)= (1/weights(i))*L_weighted.betahats.row(i);
  }
  return L_weighted;
}

partial_lasso_output partial_lasso_weighted(const arma::mat& X, const arma::colvec& y, const arma::uvec& H, const bool& partial, const arma::vec& weights, const arma::vec& grid,
                                            const double& opt_threshold, const int& opt_type){ //identical to partial_lasso(), except it calls lasso_weighted() in place of lasso().
  unsigned int T_=X.n_rows;
  unsigned int N=X.n_cols;
  unsigned int gridsize=grid.n_elem;
  unsigned int h=H.n_elem;
  arma::mat betahats(N,gridsize);
  arma::mat betahats_1(h,gridsize);
  arma::mat betahats_2(N-h,gridsize);
  arma::mat X_1=X.cols(H);
  arma::uvec minusH=linspace<arma::uvec>(0,N-1,N); minusH.shed_rows(H);
  arma::mat X_2=X.cols(minusH);
  lasso_output L;
  if(partial==true && h>0){
    arma::mat X1X1inv=inv(X_1.t()*X_1);
    arma::mat M_X1=(mat(T_,T_,fill::eye)-X_1*X1X1inv*X_1.t());
    arma::mat resX_2=M_X1*X_2;
    arma::vec resy=M_X1*y;
    L=lasso_weighted(resX_2, resy, grid, weights,
                     opt_threshold, opt_type);
    betahats_2=L.betahats;
    arma::uvec i_;
    for(unsigned int i=0; i<gridsize; i++){
      betahats_1.col(i)=X1X1inv*X_1.t()*(y-X_2*betahats_2.col(i));
      i_=i;
      betahats.submat(H,i_)=betahats_1.col(i);
      betahats.submat(minusH,i_)=betahats_2.col(i);
    }
  }else{
    L=lasso_weighted(X, y, grid, weights, opt_threshold, opt_type);
    if(h>0){
      betahats_1=L.betahats.rows(H);
      betahats_2=L.betahats.rows(minusH);
      betahats=L.betahats;
    }else{
      betahats_2=betahats=L.betahats;
    }
  }
  partial_lasso_output ret;
  ret.partial=partial;
  ret.N=N;
  ret.T_=T_;
  ret.gridsize=gridsize;
  ret.h=h;
  ret.grid=grid;
  ret.y=y;
  ret.H=H;
  ret.minusH=minusH;
  ret.betahats=betahats;
  ret.betahats_1=betahats_1;
  ret.betahats_2=betahats_2;
  ret.X=X;
  ret.X_1=X_1;
  ret.X_2=X_2;
  ret.opt_type=opt_type;
  return(ret);
} 

VAR_select_out selectTF(const arma::vec lambda_grid, const arma::mat Y, const arma::mat VAR_lags, const int& p, const int& k, 
                        const double& eps, const int& pen, const bool& pen_own, const bool& only_lag1, const double& c, 
                        const unsigned int K, double improvement_thresh, unsigned int B, const double& alpha, const double& tk){
  // c constant for plug-in approach. Paper says to take 1.1 but we can play around with it
  // K is the number of iterations of the plug-in approach. Paper says 15.
  // improvement_thresh if the % change in lambda is less than this, stop iterating. Not used in this selection method
  // B number of simulations used to estimate the quantiles of the Gaussian maximum. Not used in this selection method
  // alpha: alpha quantile of the Gaussian maximum to be considered. Overwritten in this method.
  
  
  const unsigned int& N = VAR_lags.n_cols; // N is number of predictor variables (pxd in VAR)
  const unsigned int& T_ = VAR_lags.n_rows; // number of time points
  // set up the lambda that is always used
  const double gamma_n = 0.1 / log(double(max(int(T_), p*k)));
  const double gaussian_quantile = stdnormal_inv(1-gamma_n/double(2.0*pow(k,2)*p));
  arma::vec lambda_star(1); lambda_star(0) = c*gaussian_quantile/double(sqrt(double(T_)));
  
  arma::mat out_resid(T_, k);
  arma::mat out_coef(p*k, k);
  
  arma::vec v_hat(p*k+1);
  arma::vec e2Z2(T_);
  arma::vec residual(T_);
  partial_lasso_output PLO;
  
  arma::vec constant(T_, fill::ones);
  arma::mat X=join_horiz(constant, VAR_lags); //the regressors are the same in each equation. Includes a constant
  
  for(unsigned int eq_ind=0; eq_ind<k; eq_ind++){
    // dependent variable in this equation
    arma::vec y=Y.col(eq_ind);
    // set up which variables are unpenalized in this equation
    arma::uvec H;
    if(!pen_own && only_lag1){
      arma::uvec H_temp(2); H_temp(0)=0; H_temp(1)=1+eq_ind; 
      H=H_temp;
    }else if(!pen_own && !only_lag1){
      arma::uvec H_temp(1+p); H_temp(0)=0;
      for(unsigned int lag=0; lag<p; lag++){
        H_temp(1+lag)=1+lag*k+eq_ind;
      }
      H=H_temp;
    }else{
      arma::uvec H_temp(1); H_temp(0)=0;
      H=H_temp;
    }
    residual=y; // in the initial setup, the weights are determined using the dependent variable instead of a residual
    for(unsigned int it=0; it<=K; it++){ // this loop should be done K+1 times, because it includes the initial setup
      for(unsigned int j=0; j<k*p+1; j++){
        for(unsigned int t=0; t<T_; t++){
          e2Z2(t)=pow(residual(t),2)*pow(X(t,j),2);
        }
        v_hat(j)=sqrt(mean(e2Z2));
      }
      PLO=partial_lasso_weighted(X, y, H, true, v_hat, lambda_star, eps, 3); 
      residual = y - X*PLO.betahats; 
    }
    out_resid.col(eq_ind)=residual;
    out_coef.col(eq_ind)=(PLO.betahats).submat(1,0,p*k,0);
  }
  
  VAR_select_out out;
  out.coef = out_coef;
  out.resid = out_resid;
  out.lambda = lambda_star(0);
  out.lambdas = lambda_grid;
  
  return out; 
}
// end of added by Robert

VAR_select_out selectIC(const arma::mat Y, const arma::mat VAR_lags, const int& p,  const double& nbr_lambdas,
                        const arma::vec lambda_grid, const double& eps, const int& selection, const int& pen, const bool& pen_own,
                        const bool& only_lag1, const double& tk){
  
  const int k = Y.n_cols;
  const int tt = Y.n_rows;
  
  // Non-plug-in approach
  arma::cube Phis(k, k*p, nbr_lambdas, fill::zeros);
  arma::cube resids(tt, k, nbr_lambdas, fill::zeros);
  arma::vec ics(nbr_lambdas, fill::zeros);
  arma::colvec lambda;
  
  VAR_out HVAR_fit;
  for (int il = 0; il < nbr_lambdas; il++){
    lambda = lambda_grid(il);
    HVAR_fit =  HVAR(Y, VAR_lags.t(), p, k, lambda, eps, pen, pen_own, only_lag1, tk);
    Phis.slice(il) = HVAR_fit.coef;
    resids.slice(il) = HVAR_fit.resid;
    ics(il) = VAR_ic(HVAR_fit.resid, HVAR_fit.coef, selection);
  }
  
  const int lambda_opt_index = ics.index_min();
  const arma::mat coef = Phis.slice(lambda_opt_index);
  double lambda_opt = 0;
  lambda_opt = lambda_grid[lambda_opt_index];
  
  VAR_select_out out;
  out.coef = coef.t();
  out.resid = resids.slice(lambda_opt_index);
  out.lambda = lambda_opt;
  out.lambdas = lambda_grid;
  
  return out;
}

VAR_out sparseVAR(arma::mat Y, const int& p, const bool& trim, const int& pen, const double& nbr_lambdas, //step 3 in bootstrap algorithm: estimation of the VAR by lasso
                  const double& lambda_ratio, const double& eps, const int& selection, 
                  const bool& pen_own, const bool& only_lag1,
                  const double& c, 
                  const unsigned int K, double improvement_thresh, unsigned int Nsim, const double& alpha){
  // Y: matrix (tt) x k with tt observations and k the number of variables
  // p: VAR lag order
  // trim: boolean: false fills up missing lags with zeros, thereby preserving the original time series length  for the residuals
  // pen: integer, 1 for L1 penalization; 2 for HLag penalization
  // nbr_lambdas : double, number of sparsity parameters to consider in grid (for simplicity set as double)
  // lambda_ratio : double, ratio lambda_max/lambda_min
  // eps: double, convergence tolerance
  // selection: integer, 1 for bic, 2 for aic, 3 for hq and 4 for the plug-in approach to select the sparsity parameter lambda
  // pen_own: boolean: true if own lags are to be penalized, false if they should be unpenalized
  // only_lag1 : boolean, only relevant if pen_own = false : TRUE if only first lag should be unpenalized, FALSE all own lags should be unpenalized
  // c : constant for plug-in approach (0.4 or 0.8), only relevant if selection=4
  // K : number of iterations plug-in approach, only relevant if selection=4
  // improvement_tresh: treshhold improvement plug-in approach, only relevant if selection=4
  // Nsim: number of simulations plug-in approach, only relevant if selection=4
  // alpha: quantile for plug-in approach, only relevant if selection=4
  
  const int k = Y.n_cols;
  int tt = Y.n_rows;
  
  const arma::mat VAR_lags = lag_matrix(Y, p, trim); //step 3 in bootstrap algorithm, lags are padded with 0
  
  if(trim){
    Y = Y.submat(p, 0, tt-1, k-1);
    tt = Y.n_rows;
  }
  
  // Computation of the step size
  arma::mat ZMean = mean(VAR_lags);
  arma::mat trainZt = VAR_lags;
  arma::mat Z = zeros(tt, k*p);
  for (int i = 0; i < tt; i++) {
    Z.row(i) = VAR_lags.row(i)  - ZMean ;
  }
  Z = Z.t();
  
  // Step size
  vec eigval;
  arma::mat eigvec;
  const arma::mat Zt=Z*trans(Z);
  eig_sym(eigval, eigvec, Zt);
  double tk = 1/max(eigval);
  // End computation of the step size

  
  const arma::vec lambda_grid = LambdaGridE(lambda_ratio, nbr_lambdas, Y, VAR_lags.t(), pen, p, k, TRUE, TRUE, tk); // always penalize own lags to determine lambda grid
  
  VAR_select_out VARselection;
  
  if(selection==4){// Plug-in approach
    VARselection = selectPI(lambda_grid, Y, VAR_lags, p, k, eps, pen, pen_own, only_lag1, c, K, improvement_thresh, Nsim, alpha, tk); //step 4 in bootstrap algorithm: 
  }else if(selection==5){// Theoretically founded approach
    VARselection = selectTF(lambda_grid, Y, VAR_lags, p, k, eps, pen, pen_own, only_lag1, c, K, improvement_thresh, Nsim, alpha, tk);
  }else{ // IC approach
    VARselection = selectIC(Y, VAR_lags, p, nbr_lambdas, lambda_grid, eps, selection, pen, pen_own, only_lag1, tk);
  }
  
  VAR_out out; 
  out.coef = VARselection.coef;
  out.resid = VARselection.resid;
  ///////////////////////remove
  out.lambda = VARselection.lambda;
  out.lambdas = VARselection.lambdas;
  ///////////////////////////
  return out;
}

// [[Rcpp::export]]
Rcpp::List sparseVAR_R(arma::mat Y, const int& p, const bool& trim, const int& pen, const double& nbr_lambdas,
                       const double& lambda_ratio, const double& eps, const int& selection, const double& c=0.8, 
                       const unsigned int K = 15, double improvement_thresh = 0.01, unsigned int Nsim=1000, const double& alpha = 0.05,
                       const bool& pen_own = true, const bool& only_lag1 = false){
  // Y: matrix (tt) x k with tt observations and k the number of variables
  // p: VAR lag order
  // trim: boolean: false fills up missing lags with zeros, thereby preserving the original time series length  for the residuals
  // pen: integer, 1 for L1 penalization; 2 for HLag penalization
  // nbr_lambdas : double, number of sparsity parameters to consider in grid (for simplicity set as double)
  // lambda_ratio : double, ratio lambda_max/lambda_min
  // eps: double, convergence tolerance
  // selection: integer, 1 for bic, 2 for aic, 3 for hq and 4 for the plug-in approach to select the sparsity parameter lambda
  // c : constant for plug-in approach (0.4 or 0.8), only relevant if selection=4
  // K : number of iterations plug-in approach, only relevant if selection=4
  // improvement_tresh: treshhold improvement plug-in approach, only relevant if selection=4
  // Nsim: number of simulations plug-in approach, only relevant if selection=4
  // alpha: quantile for plug-in approach, only relevant if selection=4
  // pen_own: boolean: true if own lags are to be penalized, false if they should be unpenalized
  // only_lag1 : boolean, only relevant if pen_own = false : TRUE if only first lag should be unpenalized, FALSE all own lags should be unpenalized
  
  VAR_out out = sparseVAR(Y, p, trim, pen, nbr_lambdas, lambda_ratio, eps, selection, pen_own, only_lag1, c, K, improvement_thresh, Nsim, alpha);
  
  return Rcpp::List::create(
    Rcpp::Named("coef") = out.coef,
    Rcpp::Named("resid") = out.resid
    //////////////////////////remove
    ,Rcpp::Named("lambda") = out.lambda,
    Rcpp::Named("lambdas") = out.lambdas
    ////////////////////////////////
  );
}
