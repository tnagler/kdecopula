#include <RcppArmadillo.h>
#include <kernels.h>
using namespace Rcpp;
using namespace arma;

// EPANECHNIKOV KERNELS ---------------------------------------

// [[Rcpp::export]]
arma::vec kern_epan_1d(const arma::vec& x, const double& b)
{
    vec out(x);
    int n = x.n_rows;
    
    for (int i = 0; i < n; ++i) {
        if (std::fabs(x[i] / b) >= 1.0) {
            out[i] = 0;
        } else {
            out[i] = 3.0/4.0 * (1 - pow(x[i] / b, 2));
        }
    }
    return out / b;
}

// [[Rcpp::export]]
arma::vec kern_epan_2d(const arma::vec& x, const arma::vec& y, const double& b)
{
    return kern_epan_1d(x, b) % kern_epan_1d(y, b);
}

// [[Rcpp::export]]
arma::vec kern_epan(const arma::mat& x, const arma::vec& b)
{
    vec out(x.n_rows);
    int d = x.n_cols;
    
    out = kern_epan_1d(x.col(0), b[0]);
    for (int j = 1; j < d; ++j) {
        out = out % kern_epan_1d(x.col(j), b[j]);
    }
    return out;
}

// GAUSSIAN KERNELS ------------------------------------------

// [[Rcpp::export]]
arma::vec kern_gauss_1d(const arma::vec& x, const double& b)
{
    vec out(x);
    int n = x.n_rows;
    
    for (int i = 0; i < n; ++i) {
        if (std::fabs(x[i] / b) >= 5) {
            out[i] = 0;
        } else {
            out[i] = exp(- 0.5 * pow(x[i] / b, 2)) / (sqrt(2 * M_PI)) / 0.9999994267 / b;
        }
    }
    return out;
}

// [[Rcpp::export]]
arma::vec kern_gauss_2d(const arma::vec& x, const arma::vec& y, const double& b)
{
    return kern_gauss_1d(x, b) % kern_gauss_1d(y, b);
}

// [[Rcpp::export]]
arma::vec kern_gauss(const arma::mat& x, const arma::vec& b)
{
    vec out(x.n_rows);
    int d = x.n_cols;
    
    out = kern_gauss_1d(x.col(0), b[0]);
    for (int j = 1; j < d; ++j) {
        out = out % kern_gauss_1d(x.col(j), b[j]);
    }
    return out;
}
