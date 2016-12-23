#ifndef _kdecopula_INTEGRATE_H
#define _kdecopula_INTEGRATE_H

#include <RcppArmadillo.h>
#include <interp.h>

double int_on_grid(const double& upr, 
                   const Rcpp::NumericVector& vals,
                   const Rcpp::NumericVector& grid); \

double inv_int_on_grid(const double& qq, 
                       const Rcpp::NumericVector& vals,
                       const Rcpp::NumericVector& grid);
#endif
