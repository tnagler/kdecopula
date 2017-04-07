#ifndef _kdecopula_EVALUATORS_H
#define _kdecopula_EVALUATORS_H

#include <RcppArmadillo.h>

arma::vec eval_mr(const arma::mat& uev, const arma::mat& dat, const double& b); 
arma::vec eval_beta(const arma::mat& uev, const arma::mat& dat, const double& b);
arma::vec eval_trafo(const arma::mat& uev, const arma::mat& dat, const arma::mat& B);

#endif
