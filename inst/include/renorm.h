#ifndef _kdecopula_RENORM_H
#define _kdecopula_RENORM_H

#include <RcppArmadillo.h>

Rcpp::NumericVector renorm(const Rcpp::NumericVector& x,
                           const Rcpp::NumericVector& grid, 
                           const int times,
                           const Rcpp::IntegerMatrix& helpind) ;

Rcpp::NumericVector ren_subs(const Rcpp::NumericVector& vals, 
                             const Rcpp::NumericVector& grid, 
                             const int& d); 
#endif
