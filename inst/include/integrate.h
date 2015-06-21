#ifndef _kdecopula_INTEGRATE_H
#define _kdecopula_INTEGRATE_H

#include <RcppArmadillo.h>
#include <interp.h>

double int_on_grid(double upr, Rcpp::NumericVector vals, Rcpp::NumericVector grid);
double inv_int_on_grid(double q, Rcpp::NumericVector vals, Rcpp::NumericVector grid);

#endif
