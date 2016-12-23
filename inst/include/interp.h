#ifndef _kdecopula_INTERP_H
#define _kdecopula_INTERP_H

#include <RcppArmadillo.h>

double cubic_poly(const double& x, 
                  const Rcpp::NumericVector& a); 

double cubic_indef_integral(const double& x,
                            const Rcpp::NumericVector& a); 

double cubic_integral(const double& lowr,
                      const double& upr,
                      const Rcpp::NumericVector& a);

double inv_cubic_integral(const double& q, 
                          const Rcpp::NumericVector& a);

Rcpp::NumericVector coef(const Rcpp::NumericVector& vals,
                         const Rcpp::NumericVector& grid,
                         Rcpp::NumericVector a); 

double interp_on_grid(const double& x,
                      const Rcpp::NumericVector& vals,
                      const Rcpp::NumericVector& grid,
                      Rcpp::NumericVector a);

Rcpp::NumericVector interp_2d(const Rcpp::NumericMatrix& x,
                              const Rcpp::NumericMatrix& vals, 
                              const Rcpp::NumericVector& grid,
                              Rcpp::NumericVector tmpgrid,
                              Rcpp::NumericVector tmpvals);

Rcpp::IntegerVector get(const Rcpp::IntegerMatrix& ind, 
                        const Rcpp::IntegerVector& dims); 

Rcpp::NumericVector interp(const Rcpp::NumericMatrix& x,
                           const Rcpp::NumericVector& vals, 
                           const Rcpp::NumericVector& grid, 
                           const Rcpp::IntegerMatrix& helpind); 

#endif
