#define ARMA_DONT_USE_OPENMP 1
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppParallel.h>
#include <RcppThread.h>
#include "sparseVAR.h"

using namespace arma;

struct progress {
private:
  const int max;
  const bool show_progress;
  int counter;
  int step_counter;
  std::thread::id main_id;
  tthread::mutex m;
  arma::uvec steps = arma::linspace<arma::uvec>(0, max, 21);
  
public:
  progress(const int max, const bool show_progress) : max(max), show_progress(show_progress),
  counter(0), step_counter(0), main_id(std::this_thread::get_id())
  {
    if (show_progress) {
      RcppThread::Rcout << "Progress: |------------------| \n";
      RcppThread::Rcout << "          ";
    }
  };
  
  void increment() {
    tthread::lock_guard<tthread::mutex> guard(m);
    counter++;
    if (show_progress) {
      if (std::this_thread::get_id() == main_id) {
        RcppThread::checkUserInterrupt();
        if (counter > steps(step_counter + 1)) {
          RcppThread::Rcout << "*";
          step_counter++;
        }
      }
    }
  }
  
  ~progress() {
    if (show_progress) {
      RcppThread::Rcout << "*\n";
    }
  }
};

struct boot_out {
  double mean;
  arma::vec boot_quantiles;
  int par;
  //////remove 
  arma::mat coef_pre;
  arma::mat coef_post;
  arma::vec lambdas;
  double lambda;
  arma::mat smeans;
  ////////////
};

arma::mat custom_rnorm(const unsigned int& r, const unsigned int& c,
                       const double& mu = 0, const double& sd = 1) {
  Rcpp::NumericVector draw = Rcpp::rnorm(r * c, mu, sd);
  const arma::mat Z = arma::mat(draw.begin(), r, c, false, true);
  // const arma::mat Z = arma::randn(r, c, distr_param(mu, sd));
  return (Z);
}

arma::mat custom_runif(const unsigned int& r, const unsigned int& c,
                       const double& min = 0, const double& max = 1) {
  Rcpp::NumericVector draw = Rcpp::runif(r * c);
  const arma::mat U = arma::mat(draw.begin(), r, c, false, true);
  // const arma::mat U = arma::randu(r, c, distr_param(min, max));
  return (U);
}

arma::umat custom_sample(const unsigned int& r, const unsigned int& c,
                         const unsigned int& n) {
  // Rcpp::NumericVector draw = Rcpp::sample(n, r * c, true) - 1;
  // const arma::umat i = arma::mat(draw.begin(), r, c, false, true);
  // const arma::umat i = arma::randi(r, c, distr_param(0, n - 1));
  const arma::uvec i_vec = Rcpp::RcppArmadillo::sample(linspace<arma::uvec>(0, n - 1, n), r * c, true);
  const arma::umat i = reshape(i_vec, r, c);
  return (i);
}

arma::mat diff(const arma::mat& x, const bool& trim = true, const double& c = 1.0) {
  const arma::mat d_x = x - c * lag_matrix(x, 1, false);
  return d_x.rows(trim, x.n_rows - 1);
}

arma::mat ols(const arma::mat& y, const arma::mat& x) {
  const arma::mat b = arma::inv_sympd(x.t() * x) * x.t() * y;
  return b;
}

VAR_out VAR(const arma::mat& y, const int& p, const bool& intercept = false) {
  arma::mat lags_y = lag_matrix(y, p, false);
  // Setting trim false fills up missing lags with zeros, thereby preserving
  // the original time series length for the residuals
  VAR_out out;
  if (intercept) {
    lags_y = join_horiz(ones(lags_y.n_rows), lags_y);
  }
  arma::mat b = ols(y, lags_y);
  out.coef = b.tail_rows(lags_y.n_cols - intercept);
  out.resid = y - lags_y * b;
  return out;
}

arma::sp_mat companion_form(const arma::mat A) {
  const unsigned int k = A.n_cols;
  const unsigned int p = A.n_rows / k;
  arma::sp_mat B(p * k, p * k);
  B.rows(0, k - 1) = A.t();
  if (p > 1) {
    B.submat(k, 0, k * p - 1, (p - 1) * k - 1) = speye((p - 1) * k, (p - 1) * k);
  }
  return B;
}

void VAR_root_bound(VAR_out& V, const double& max_EV = 0.999) {
  const arma::sp_mat A = companion_form(V.coef);
  eigs_opts opts;
  //  opts.maxiter = 10000;
  //  opts.subdim = 50;
  //  opts.tol = 0.0000001;
  //  arma::cx_vec eigval = eigs_gen(A, 1, "lm", opts);
  arma::cx_vec eigval = eig_gen(arma::mat(A));
  const arma::vec A_EV = abs(eigval);
  if (A_EV(0) > max_EV) {
    V.coef = max_EV * V.coef / A_EV(0);
  }
}

// [[Rcpp::export]]
int VAR_determine_p(const arma::mat& x, const int& pmax = 10, const int& criterion = 1){
  // criterion : 1 bic; 2 aic; 3 hq
  
  const int n = x.n_rows;
  const int k = x.n_cols;
  double C_T; 
  
  if (criterion == 1){ //BIC
    C_T = log(n); 
  } else if (criterion == 2){ //AIC
    C_T = 2;
  } else if (criterion == 3){ //HQ
    C_T = 2*log(log(n));
  }
  
  VAR_out fit; 
  arma::mat resids(n, k, fill::zeros);
  arma::mat omega; 
  arma::vec IC(pmax, fill::zeros);
  
  for (int ip = 0; ip < pmax; ip++) {
    
    for(int ik=0; ik < k; ik++){
      fit = VAR(x.col(ik), ip, true);
      resids.submat(0, ik, n-1, ik) = fit.resid; 
    }
    
    omega = diagmat(resids.t() * resids)/n;  
    IC(ip) = log(det(omega)) + C_T*ip*k/n; 
  }
  
  const int popt = IC.index_min() + 1;
  
  return popt;
}

// [[Rcpp::export]]
double determine_block_length(const arma::mat& x){
  const int n = x.n_rows;
  const int k = x.n_cols;
  VAR_out fit; 
  arma::vec rho(k, fill::zeros);
  arma::vec sigma(k, fill::zeros);
  
  for(int ik = 0; ik < k; ik++) {
    fit = VAR(x.col(ik), 1, true);
    sigma(ik) = var(vec(fit.resid));
    rho(ik) = fit.coef(0, 0);
  }
  arma::vec a1n = 4 * square(rho % sigma) / (pow(1 - rho, 6) % pow(1 + rho, 2));
  arma::vec a1d = square(sigma) / pow(1 - rho, 4);
  double a1 = mean(a1n) / mean(a1d);
  double S = 1.1447 * pow(a1 * n, 0.3333333);
  
  return S;
}


arma::mat sorted_means(const arma::mat& y, const bool& abs_val = true) {
  const unsigned int N = y.n_cols;
  arma::mat sort_means = arma::zeros(N, 2);
  arma::rowvec means = mean(y, 0);
  if (abs_val) {
    means = abs(means);
  }
  const arma::uvec i = sort_index(means, "descend");
  sort_means.col(0) = means.elem(i);
  sort_means.col(1) = linspace(0, N - 1, N).elem(i);
  return sort_means;
}

// [[Rcpp::export]]
arma::mat gen_VAR(const arma::mat& u, const arma::mat& ar,
                  const arma::mat& init, const bool& include_init = false){
  const int T = u.n_rows;
  const int k = ar.n_cols;
  const int p = ar.n_rows / k;
  arma::mat y = zeros(T + p, k);
  if (init.n_rows == p & init.n_cols == k) {
    y.rows(0, p - 1) = init;
  } else {
    y.rows(0, p - 1).fill(init(0, 0));
  }
  arma::rowvec lags_y;
  for (int iT = 0; iT < T; iT++) { //step 8 of the bootstrap algorithm
    lags_y = lag_matrix(y.rows(iT, iT + p), p, false).tail_rows(1);
    y.row(iT + p) = lags_y * ar + u.row(iT);
  }
  if (!include_init) {
    y = y.tail_rows(T);
  }
  return y;
}

//[[Rcpp::export]]
arma::mat sim_mvn_chol(const arma::mat& Sigma, const unsigned int& T) {
  unsigned int N = Sigma.n_rows;
  const arma::sp_mat Sigma_sqrt = sp_mat(chol(Sigma));
  const arma::mat Z = custom_rnorm(T, N);
  const arma::mat Y = Z * Sigma_sqrt;
  return Y;
}

// [[Rcpp::export]]
arma::mat sim_VAR_cpp(const unsigned int& T, const arma::mat& ar, const arma::mat& Sigma,
                      const unsigned int& burn){
  const arma::mat u = sim_mvn_chol(Sigma, T + burn);
  const arma::mat x = gen_VAR(u, ar, zeros(1, 1));
  return x.tail_rows(T);
}

// [[Rcpp::export]]
Rcpp::List sim_VAR_cpp_both(const unsigned int& T, const arma::mat& ar, const arma::mat& Sigma,
                          const unsigned int& burn){
  const arma::mat u = sim_mvn_chol(Sigma, T + burn);
  const arma::mat x = gen_VAR(u, ar, zeros(1, 1));
  return Rcpp::List::create(Rcpp::Named("x")=x.tail_rows(T),
                            Rcpp::Named("u")=u.tail_rows(T)
  );
}

arma::mat SB(const arma::mat& e, const arma::uvec& i, 
             const arma::mat& ar_est, const arma::mat& init){
  const int T = e.n_rows;
  const arma::uvec index = i.subvec(0, T - 1);
  const arma::mat e_star = e.rows(index);
  arma::mat x_star = gen_VAR(e_star, ar_est, init, false);
  return x_star.tail_rows(T);
}

arma::mat SWB(const arma::mat& e, const arma::vec& z, 
              const arma::mat& ar_est, const arma::mat& init){
  const int T = e.n_rows;
  const int k = e.n_cols;
  const arma::mat e_star = repelem(z, 1, k) % e; //step 7 of the bootstrap algorithm
  arma::mat x_star = gen_VAR(e_star, ar_est, init, false);
  return x_star.tail_rows(T);
}

arma::mat MBB(const arma::mat& x, const arma::uvec& i, const int& l){
  const int T = x.n_rows;
  const int N = x.n_cols;
  const int nb = ceil(double(T) / double(l));
  arma::mat x_star = zeros(nb*l, N);
  for(int j = 0; j < nb; j++){
    x_star.rows(j*l, j*l + l - 1) = x.rows(i(j), i(j) + l - 1);
  }
  return x_star.head_rows(T);
}

arma::mat BWB(const arma::mat& x, const arma::vec& z, const int& l){
  const int T = x.n_rows;
  const int N = x.n_cols;
  const arma::mat xi_rep = repelem(z, l, N);
  const arma::mat x_star = x % xi_rep.head_rows(T);
  return x_star;
}

struct boot_sample_VAR_SB : public RcppParallel::Worker
{
  // inputs
  const arma::mat& e;
  const arma::umat& i;
  const arma::mat& ar;
  const arma::mat& smeans;
  const bool& abs_val;
  const arma::mat& init;
  
  const unsigned int N = e.n_cols;
  const unsigned int p = ar.n_rows / N;
  const unsigned int B = i.n_cols;
  
  // Output
  arma::cube& means_boot;
  arma::cube& x_boot;
  progress& prog;
  
  // initialize with source and destination
  boot_sample_VAR_SB(const arma::mat& e, const arma::umat& i, const arma::mat& ar, 
                     const arma::mat& smeans, const bool& abs_val, const arma::mat& init, 
                     arma::cube& means_boot, arma::cube& x_boot, progress& prog)
    : e(e), i(i), ar(ar), smeans(smeans), abs_val(abs_val),
      init(init), means_boot(means_boot), x_boot(x_boot), prog(prog) {}
  
  // Bootstrap
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t iB = begin; iB < end; iB++) {
      x_boot.slice(iB) = SB(e, i.col(iB), ar, init);
      means_boot.slice(iB) = sorted_means(x_boot.slice(iB), abs_val);
      prog.increment();
    }
  }
};

struct boot_sample_VAR_SWB : public RcppParallel::Worker
{
  // inputs
  const arma::mat& e;
  const arma::mat& z, ar;
  const arma::mat& smeans;
  const bool& abs_val;
  const arma::mat& init;
  
  const unsigned int N = e.n_cols;
  const unsigned int p = ar.n_rows / N;
  const unsigned int B = z.n_cols;
  
  // Output
  arma::cube& means_boot;
  arma::cube& x_boot;
  progress& prog;
  
  // initialize with source and destination
  boot_sample_VAR_SWB(const arma::mat& e, const arma::mat& z, const arma::mat& ar, 
                      const arma::mat& smeans, const bool& abs_val, const arma::mat& init, 
                      arma::cube& means_boot, arma::cube& x_boot, progress& prog)
    : e(e), z(z), ar(ar), smeans(smeans), abs_val(abs_val),
      init(init), means_boot(means_boot), x_boot(x_boot), prog(prog) {}
  
  // Bootstrap
  void operator()(std::size_t begin, std::size_t end) { 
    for (std::size_t iB = begin; iB < end; iB++) { //step 5 of the bootstrap algorithm
      x_boot.slice(iB) = SWB(e, z.col(iB), ar, init);
      means_boot.slice(iB) = sorted_means(x_boot.slice(iB), abs_val);
      prog.increment();
    }
  }
};

struct boot_sample_MBB : public RcppParallel::Worker
{
  // inputs
  const arma::mat& x;
  const arma::umat& i;
  const int l;
  const arma::mat& smeans;
  const bool& abs_val;
  
  const unsigned int N = x.n_cols;
  const unsigned int B = i.n_cols;
  
  // Output
  arma::cube& means_boot;
  arma::cube& x_boot;
  progress& prog;
  
  // initialize with source and destination
  boot_sample_MBB(const arma::mat& x, const arma::umat& i, 
                  const int& l, const arma::mat& smeans, const bool& abs_val,
                  arma::cube& means_boot, arma::cube& x_boot, progress& prog)
    : x(x), i(i), l(l), smeans(smeans), abs_val(abs_val),
      means_boot(means_boot), x_boot(x_boot), prog(prog) {}
  
  // Bootstrap
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t iB = begin; iB < end; iB++) {
      x_boot.slice(iB) = MBB(x, i.col(iB), l);
      means_boot.slice(iB) = sorted_means(x_boot.slice(iB), abs_val);
      prog.increment();
    }
  }
};

struct boot_sample_BWB : public RcppParallel::Worker
{
  // inputs
  const arma::mat& x;
  const arma::mat& z;
  const int l;
  const arma::mat& smeans;
  const bool& abs_val;
  
  const unsigned int N = x.n_cols;
  const unsigned int B = z.n_cols;
  
  // Output
  arma::cube& means_boot;
  arma::cube& x_boot;
  progress& prog;
  
  // initialize with source and destination
  boot_sample_BWB(const arma::mat& x, const arma::mat& z,
                  const int& l, const arma::mat& smeans, const bool& abs_val,
                  arma::cube& means_boot, arma::cube& x_boot, progress& prog)
    : x(x), z(z), l(l), smeans(smeans), abs_val(abs_val),
      means_boot(means_boot), x_boot(x_boot), prog(prog) {}
  
  // Bootstrap
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t iB = begin; iB < end; iB++) {
      x_boot.slice(iB) = BWB(x, z.col(iB), l);
      means_boot.slice(iB) = sorted_means(x_boot.slice(iB), abs_val);
      prog.increment();
    }
  }
};

boot_out boot_means(const arma::mat& x_with_means, const arma::mat oracle_A, const arma::mat oracle_u, const int& boot, const int& p, const int& l, 
                    const bool& abs_val, const arma::vec& q, const int& B, const arma::mat init, 
                    const bool& show_progress, const int& penalization, 
                    const double& nbr_lambdas, const double& lambda_ratio, 
                    const int& selection, const double& eps, const bool& pen_own, const bool& only_lag1,
                    const double& c, 
                    const unsigned int K, double improvement_thresh, unsigned int Nsim, const double& alpha) {
  // x_with_means: this is the raw data, once demeaned we simply call it x
  // penalization: integer, 0 (no penalization), 1 (L1), 2 (HLag)
  // nbr_lambdas : double, number of sparsity parameters to consider in grid (for simplicity set as double)
  // lambda_ratio : double, ratio lambda_max/lambda_min
  // selection: integer, 1 for bic, 2 for aic, 3 for hq to select the sparsity parameter lambda
  // eps: double, convergence tolerance
  // pen_own: boolean: true if own lags are to be penalized, false if they should be unpenalized
  // only_lag1 : boolean, only relevant if pen_own = false : TRUE if only first lag should be unpenalized, FALSE all own lags should be unpenalized
  
  
  arma::mat x=x_with_means;
  int T = x.n_rows;
  int N = x.n_cols;
  
  //demean x: step 2 in the bootstrap algorithm
  double m=0.0; 
  for(unsigned int i=0; i<N; i++){
    m = mean(x_with_means.col(i));
    x.col(i)=x_with_means.col(i)-m;
  }
  
  VAR_out out;
  int p_boot, l_boot, nb;
  boot_out out_boot;
  
  if (boot == 1 | boot == 2){
    if (p == 0) {
      p_boot = VAR_determine_p(x);
    } else {
      p_boot = p;
    }
    out_boot.par = p_boot;
    ///add a case here where the out is built from the true thing
    if (penalization == 0){
      out = VAR(x, p_boot);
    } else if (penalization == 1){
      out = sparseVAR(x, p_boot, false, 1, nbr_lambdas, lambda_ratio, eps, selection, pen_own, only_lag1, c, K, improvement_thresh, Nsim, alpha); //step 3 in bootstrap algorithm
    } else if (penalization == 2){
      out = sparseVAR(x, p_boot, false, 2, nbr_lambdas, lambda_ratio, eps, selection, pen_own, only_lag1, c, K, improvement_thresh, Nsim, alpha);
    } else if (penalization == -1){ // oracle
      out.coef=oracle_A; 
      out.resid=oracle_u;
    }
    
    ////////////////////remove
    out_boot.coef_pre=out.coef;
    //////////////////////
    
    VAR_root_bound(out);
    
    ////////////////////remove
    out_boot.coef_post=out.coef;
    //////////////////////
  } else {
    if (l == 0) {
      l_boot = round(determine_block_length(x));
    } else {
      l_boot = l;
    }
    l_boot = std::max(1, std::min(l_boot, T/2));
    out_boot.par = l_boot;
  }
  
  arma::cube means_boot(N, 2, B);
  arma::cube x_boot(T, N, B);
  
  // this part needs to look at the original data to compute the statistic
  arma::mat smeans = sorted_means(x_with_means, abs_val); //step 1 in the bootstrap algorithm
  ////////////////////remove
  out_boot.smeans=smeans;
  //////////////////////
 
  
  
  progress prog(B, show_progress);
  if (boot == 1) {//
    const arma::mat z = custom_rnorm(T, B, 0, 1); //step 6 of the bootstrap algorithm
    boot_sample_VAR_SWB boot_sample_x(out.resid, z, out.coef, smeans, abs_val, init, means_boot, x_boot, prog);
    RcppParallel::parallelFor(0, B, boot_sample_x);
  } else if (boot == 2){
    const arma::umat i = custom_sample(T, B, T);
    boot_sample_VAR_SB boot_sample_x(out.resid, i, out.coef, smeans, abs_val, init, means_boot, x_boot, prog);
    RcppParallel::parallelFor(0, B, boot_sample_x);
  } else if (boot == 3){
    nb = ceil(double(T) / double(l_boot));
    const arma::mat z = custom_rnorm(nb, B, 0, 1);
    boot_sample_BWB boot_sample_x(x, z, l_boot, smeans, abs_val, means_boot, x_boot, prog);
    RcppParallel::parallelFor(0, B, boot_sample_x);
  } else if (boot == 4){
    nb = ceil(double(T) / double(l_boot));
    const arma::umat i = custom_sample(nb, B, T - l_boot + 1);
    boot_sample_MBB boot_sample_x(x, i, l_boot, smeans, abs_val, means_boot, x_boot, prog);
    RcppParallel::parallelFor(0, B, boot_sample_x);
  }
  
  arma::vec max_means_q(2);
  arma::vec max_boot_means = means_boot.subcube(0, 0, 0, 0, 0, B - 1); //step 9 in the bootstrap algorithm
  arma::vec qu = quantile(max_boot_means, q);
  
  out_boot.mean = smeans(0); //step 1 in the bootstrap algorithm
  out_boot.boot_quantiles = qu;
  return out_boot;
}

// [[Rcpp::export]]
Rcpp::List boot_means_R(const arma::mat& x, const arma::mat oracle_A, const arma::mat oracle_u, const int& boot = 1, const int& p = 1, 
                        const int& l = 1, const bool& abs_val = true, 
                        const arma::vec& q = 0.95 * ones(1), const int& B = 9999, 
                        const bool& show_progress = false, const int& penalization = 1, 
                        const double& nbr_lambdas = 10, const double& lambda_ratio = 100,
                        const int& selection = 1, const double& eps = 0.001, const bool& pen_own = true, const bool& only_lag1 = false,
                        const double& c=0.8, 
                        const unsigned int K = 15, double improvement_thresh = 0.01, unsigned int Nsim=1000, const double& alpha = 0.05) {
  boot_out out = boot_means(x, oracle_A, oracle_u, boot, p, l, abs_val, q, B, zeros(1, 1), show_progress,
                            penalization, nbr_lambdas, lambda_ratio, selection, eps, pen_own, only_lag1, c, K, improvement_thresh, Nsim, alpha);
  return Rcpp::List::create(
    Rcpp::Named("mean") = out.mean,
    Rcpp::Named("boot_quantiles") = out.boot_quantiles,
    ////////////////////////////////////////remove
    Rcpp::Named("coef_pre") = out.coef_pre,
    Rcpp::Named("coef_post") = out.coef_post,
    Rcpp::Named("lambdas") = out.lambdas,
    Rcpp::Named("lambda") = out.lambda,
    Rcpp::Named("smeans") = out.smeans,
    //////////////////////////////////////////////
    Rcpp::Named("par") = out.par);
}

// [[Rcpp::export]]
Rcpp::List VAR_R(const arma::mat& y, const int& p, const bool& intercept = true) {
  VAR_out out = VAR(y, p);
  return Rcpp::List::create(
    Rcpp::Named("coef") = out.coef,
    Rcpp::Named("resid") = out.resid);
}