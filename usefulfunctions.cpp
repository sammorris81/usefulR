// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
Rcpp::NumericVector rowSumsC(arma::mat X) {
  int nrow = X.n_rows;

  Rcpp::NumericVector out(nrow);
  for (uword i = 0; i < nrow; i++) {
    out[i] = accu(X.row(i));
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector colSumsC(arma::mat X) {
  int ncol = X.n_cols;

  Rcpp::NumericVector out(ncol);
  for (uword i = 0; i < ncol; i++) {
    out[i] = accu(X.col(i));
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sweepC1plus(Rcpp::NumericMatrix X,
                                Rcpp::NumericVector y) {
  int nrow = X.nrow(); int ncol = X.ncol();

  Rcpp::NumericMatrix out(nrow, ncol);
  for (uword i = 0; i < ncol; i++) {
    out(_, i) = X(_, i) + y;
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sweepC2plus(Rcpp::NumericMatrix X,
                                Rcpp::NumericVector y) {
  int nrow = X.nrow(); int ncol = X.ncol();

  Rcpp::NumericMatrix out(nrow, ncol);
  for (uword i = 0; i < nrow; i++) {
    out(i, _) = X(i, _) + y;
  }

  return out;
}

// [[Rcpp::export]]
arma::mat sweepC1times(arma::mat X,
                       arma::colvec y) {
  int nrow = X.n_rows; int ncol = X.n_cols;

  arma::mat out(nrow, ncol);
  for (uword i = 0; i < ncol; i++) {
    out.col(i) = X.col(i) % y;
  }

  return out;
}

// [[Rcpp::export]]
arma::mat sweepC2times(arma::mat X,
                       arma::rowvec y) {
  int nrow = X.n_rows; int ncol = X.n_cols;

  arma::mat out(nrow, ncol);
  for (uword i = 0; i < nrow; i++) {
    out.row(i) = X.row(i) % y;
  }

  return out;
}

// [[Rcpp::export]]
arma::mat powC(arma::mat X, double alpha) {
  return exp(alpha * log(X));
}

