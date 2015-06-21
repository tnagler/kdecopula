#ifndef _kdecopula_HFUNC_H
#define _kdecopula_HFUNC_H

#include <RcppArmadillo.h>
#include <interp.h>
#include <integrate.h>

Rcpp::NumericVector eval_hfunc_2d(Rcpp::NumericMatrix uev, int cond_var, Rcpp::NumericMatrix vals, Rcpp::NumericVector grid);
Rcpp::NumericVector eval_hfunc(Rcpp::NumericMatrix uev, Rcpp::IntegerVector cond_var, Rcpp::IntegerVector uncond_var,
Rcpp::NumericVector vals, Rcpp::NumericVector grid, Rcpp::NumericMatrix helpgrid, Rcpp::IntegerMatrix helpind); 
Rcpp::NumericVector inv_hfunc(Rcpp::NumericMatrix uev, int cond_var, Rcpp::NumericMatrix vals, Rcpp::NumericVector grid); 
Rcpp::NumericVector eval_cdf(Rcpp::NumericMatrix uev, Rcpp::NumericVector vals, Rcpp::NumericVector grid, 
Rcpp::NumericMatrix helpgrid, Rcpp::IntegerMatrix helpind);

#endif
