#ifndef _kdecopula_INTERP_H
#define _kdecopula_INTERP_H

#include <RcppArmadillo.h>

double cubic_poly(double x, Rcpp::NumericVector a);
double cubic_indef_integral(double x, Rcpp::NumericVector a); 
double cubic_integral(double lowr, double upr, Rcpp::NumericVector a);
double inv_cubic_integral(const double q, Rcpp::NumericVector a);
Rcpp::NumericVector coef(Rcpp::NumericVector vals, Rcpp::NumericVector grid); 
double interp_on_grid(double x, Rcpp::NumericVector vals, Rcpp::NumericVector grid);
Rcpp::NumericVector interp_2d(Rcpp::NumericMatrix x, Rcpp::NumericMatrix vals, Rcpp::NumericVector grid);
Rcpp::IntegerVector get(Rcpp::IntegerMatrix ind, Rcpp::IntegerVector dims);
Rcpp::NumericVector interp(Rcpp::NumericMatrix x, Rcpp::NumericVector vals, Rcpp::NumericVector grid, Rcpp::IntegerMatrix helpind);

#endif
