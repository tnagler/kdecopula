#ifndef _kdecopula_KERNELS_H
#define _kdecopula_KERNELS_H

#include <RcppArmadillo.h>

arma::vec kern_epan_1d(const arma::vec& x, const double& b);
arma::vec kern_epan_2d(const arma::vec& x, const arma::vec& y, const double& b);
arma::vec kern_epan(const arma::mat& x, const arma::vec& b);
arma::vec kern_gauss_1d(const arma::vec& x, const double& b);
arma::vec kern_gauss_2d(const arma::vec& x, const arma::vec& y, const double& b);
arma::vec kern_gauss(const arma::mat& x, const arma::vec& b);
    
#endif
