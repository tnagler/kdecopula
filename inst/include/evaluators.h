#ifndef _kdecopula_EVALUATORS_H
#define _kdecopula_EVALUATORS_H

#include <RcppArmadillo.h>
#include <kernels.h>

Rcpp::NumericVector eval_mr(Rcpp::NumericMatrix uev, Rcpp::NumericMatrix dat, double b);
Rcpp::NumericVector eval_beta(Rcpp::NumericMatrix uev, Rcpp::NumericMatrix dat, double b);
Rcpp::NumericVector eval_trafo(Rcpp::NumericMatrix uev, Rcpp::NumericMatrix dat, arma::mat B);

#endif
