#include <RcppArmadillo.h>
#include <kernels.h>
using namespace Rcpp;
using namespace arma;

//' Epanechnikov kernel (univariate)
//' 
//' @param x vector of evaluation points.
//' @param b bandwidth parameter.
//' 
//' @noRd
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

//' Epanechnikov product kernel (bivariate)
//' 
//' Computes \eqn{K(x/b, y/b) / b}.
//' 
//' @param x vector of evaluation points.
//' @param y vector of evaluation points.
//' @param b bandwidth parameter.
//' 
//' @noRd
// [[Rcpp::export]]
arma::vec kern_epan_2d(const arma::vec& x, const arma::vec& y, const double& b)
{
    return kern_epan_1d(x, b) % kern_epan_1d(y, b);
}

//' Epanechnikov product kernel (arbitray dimension)
//' 
//' @param x mxd matrix of evaluation points.
//' @param b length d vector of bandwidth parameters.
//' 
//' @noRd
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

//' Gaussian kernel (univariate)
//' 
//' @param x vector of evaluation points.
//' @param b bandwidth parameter.
//' 
//' @noRd
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

//' Gaussian product kernel (bivariate)
//' 
//' Computes \eqn{K(x/b, y/b) / b}.
//' 
//' @param x vector of evaluation points.
//' @param y vector of evaluation points.
//' @param b bandwidth parameter.
//' 
//' @noRd
// [[Rcpp::export]]
arma::vec kern_gauss_2d(const arma::vec& x, const arma::vec& y, const double& b)
{
    return kern_gauss_1d(x, b) % kern_gauss_1d(y, b);
}

//' Gaussian product kernel (arbitray dimension)
//' 
//' @param x mxd matrix of evaluation points.
//' @param b length d vector of bandwidth parameters.
//' 
//' @noRd
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
