#ifndef _kdecopula_RENORM_H
#define _kdecopula_RENORM_H

#include <RcppArmadillo.h>

Rcpp::NumericVector renorm(Rcpp::NumericVector x, Rcpp::NumericVector grid, int times, Rcpp::IntegerMatrix helpind);
Rcpp::NumericVector ren_subs(Rcpp::NumericVector vals, Rcpp::NumericVector grid);

#endif
