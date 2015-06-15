#ifndef _kdecopula_KERNELS_H
#define _kdecopula_KERNELS_H

#include <RcppArmadillo.h>

Rcpp::NumericVector kern_epan(Rcpp::NumericVector x);
Rcpp::NumericVector kern_epan_2d(Rcpp::NumericVector x, Rcpp::NumericVector y);
Rcpp::NumericVector kern_gauss(Rcpp::NumericVector x);
Rcpp::NumericVector kern_gauss_2d(Rcpp::NumericVector x, Rcpp::NumericVector y);

#endif
