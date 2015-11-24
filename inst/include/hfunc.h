#ifndef _kdecopula_HFUNC_H
#define _kdecopula_HFUNC_H

#include <RcppArmadillo.h>
#include <interp.h>
#include <integrate.h>

Rcpp::NumericVector eval_hfunc_2d(const Rcpp::NumericMatrix& uev,
                                  const int& cond_var, 
                                  const Rcpp::NumericMatrix& vals,
                                  const Rcpp::NumericVector& grid);

Rcpp::NumericVector eval_hfunc(const Rcpp::NumericMatrix& uev, 
                               const Rcpp::IntegerVector& cond_var,
                               const Rcpp::IntegerVector& uncond_var,
                               const Rcpp::NumericVector& vals,
                               const Rcpp::NumericVector& grid, 
                               const Rcpp::NumericMatrix& helpgrid,
                               const Rcpp::IntegerMatrix& helpind);

Rcpp::NumericVector inv_hfunc(const Rcpp::NumericMatrix& uev,
                              const int& cond_var,
                              const Rcpp::NumericMatrix& vals,
                              const Rcpp::NumericVector& grid);

Rcpp::NumericVector eval_cdf(const Rcpp::NumericMatrix& uev, 
                             const Rcpp::NumericVector& vals,
                             const Rcpp::NumericVector& grid, 
                             const Rcpp::NumericMatrix& helpgrid, 
                             const Rcpp::IntegerMatrix& helpind);

#endif
