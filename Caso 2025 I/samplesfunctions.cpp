// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec sample_theta_cpp(
    const arma::vec& k2jk,
    const arma::ivec& njk,
    double sig_beta,
    double sig_E,
    double sig_M,
    double sig_D,
    const arma::vec& beta_0,
    const arma::mat& X,
    const arma::mat& tX,
    const arma::vec& y,
    int e,
    int l,
    int d
) {
  int p = beta_0.n_elem;
  arma::vec sigma_diag(p);
  sigma_diag(0) = sig_beta;
  sigma_diag.subvec(1, e) = arma::vec(e, fill::value(sig_E));
  sigma_diag.subvec(e + 1, e + l) = arma::vec(l, fill::value(sig_M));
  sigma_diag.subvec(e + l + 1, e + l + d) = arma::vec(d, fill::value(sig_D));  
  arma::mat sigma_inv = diagmat(1.0 / sigma_diag);
  
  arma::uword N = y.n_elem;
  arma::vec sigma_vec(N);
  for (int j = 0, idx = 0; j < njk.n_elem; ++j) {
    for (int k = 0; k < njk[j]; ++k) {
      sigma_vec[idx++] = 1.0 / k2jk[j];
    }
  }
  
  arma::mat XtSigma = tX.each_row() % sigma_vec.t();
  arma::mat V = arma::inv_sympd(sigma_inv + XtSigma * X);
  arma::vec m = V * (sigma_inv * beta_0 + XtSigma * y);
  
  return mvnrnd(m, V);
}

// [[Rcpp::export]]
arma::vec sample_k2jk_cpp(const arma::vec& y,
                          double nu_k,
                          const arma::vec& nk,
                          const arma::ivec& njk,
                          int nj,
                          const arma::rowvec& theta,
                          const arma::mat& X,
                          const arma::ivec& mpios_id,
                          const arma::vec& k2_k) {
  
  int n = y.n_elem;
  arma::vec residuos = y - X * theta.t();
  arma::vec residuos2 = arma::square(residuos);
  
  arma::vec sum_res(nj, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    int j = mpios_id(i); // municipio
    sum_res(j) += residuos2(i);
  }
  
  arma::vec a = (nu_k + arma::conv_to<arma::vec>::from(njk)) / 2.0;
  arma::vec b(nj);
  for (int j = 0; j < nj; ++j) {
    int k = nk(j);  // aquí nk debe mapear municipio j → departamento k
    b(j) = (nu_k * k2_k(k) + sum_res(j)) / 2.0;
  }
  
  arma::vec k2jk(nj);
  for (int j = 0; j < nj; ++j) {
    k2jk(j) = 1.0 / R::rgamma(a(j), 1.0 / b(j));
  }
  
  return k2jk;
}