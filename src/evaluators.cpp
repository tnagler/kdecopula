#include <RcppArmadillo.h>
#include <kernels.h>
#include <interp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// MIRROR-REFLECTION ETIMATOR
// [[Rcpp::export]]
arma::vec eval_mr(const arma::mat& uev,  /* evaluation points */
                  const arma::mat& dat,  /* cpoula data       */
                  const double& b)       /* bandwidth         */
{
    double n = dat.n_rows;
    int m = uev.n_rows;
    vec out(uev.n_rows);
    vec tmp1(n), tmp2(n), tmp3(n),
    tmp4(n), tmp5(n), tmp6(n),
    tmp7(n), tmp8(n), tmp9(n);

    for (int i = 0; i < m; ++i) {
        // compute kernels at all reflections
        tmp1 = kern_epan_2d(uev(i, 0) - dat.col(0), uev(i, 1) - dat.col(1), b);
        tmp2 = kern_epan_2d(uev(i, 0) + dat.col(0), uev(i, 1) - dat.col(1), b);
        tmp3 = kern_epan_2d(uev(i, 0) - dat.col(0), uev(i, 1) + dat.col(1), b);
        tmp4 = kern_epan_2d(uev(i, 0) + dat.col(0), uev(i, 1) + dat.col(1), b);
        tmp5 = kern_epan_2d(uev(i, 0) - dat.col(0), uev(i, 1) + dat.col(1) - 2, b);
        tmp6 = kern_epan_2d(uev(i, 0) + dat.col(0), uev(i, 1) + dat.col(1) - 2, b);
        tmp7 = kern_epan_2d(uev(i, 0) + dat.col(0) - 2, uev(i, 1) - dat.col(1), b);
        tmp8 = kern_epan_2d(uev(i, 0) + dat.col(0) - 2, uev(i, 1) + dat.col(1), b);
        tmp9 = kern_epan_2d(uev(i, 0) + dat.col(0) - 2, uev(i, 1) + dat.col(1) - 2, b);
        
        // estimate is sum of all kernels divided by nb^2
        out[i] = sum(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8) / n;
        
    }
    return out;
}

// BETA KERNEL ESTIMATOR
// [[Rcpp::export]]
arma::vec eval_beta(const arma::mat& uev,  /* evaluation points */
                    const arma::mat& dat,  /* cpoula data       */
                    double b)              /* bandwidth         */
{
    int n = dat.n_rows;
    int d = dat.n_cols;
    int m = uev.n_rows;
    vec datj(n), B0(n), B1(n), out(uev.n_rows);
    
    for (int i = 0; i < m; ++i) {
        // initialize with vector of ones
        vec B = rep(1.0, n);
        
        // multiply with a beta kernel for each dimension
        for (int j = 0; j < d; ++j) {
            datj = dat.col(j);
            B0 = dbeta(as<NumericVector>(wrap(datj)), 
                       uev(i, j) / b + 1,
                       (1 - uev(i, j)) / b + 1);
            B = B % B0;
        }
        
        // copula density estimate is mean of product kernels
        out[i] = mean(B);
    }
    return out;
}


// TRANSFORMATION ESTIMATOR
// [[Rcpp::export]]
arma::vec eval_trafo(const arma::mat& uev,  /* evaluation points */
                     const arma::mat& dat,  /* cpoula data       */
                     const arma::mat& B)    /* bandwidth matrix  */
{ 
    int n = dat.n_rows;
    int d = dat.n_cols;
    int m = uev.n_rows;
    vec out(m);
    
    // transform data by inverse Gaussian cdf
    mat xev = as<vec>(wrap(qnorm(as<NumericVector>(wrap(uev))))); 
    mat xdat = as<vec>(wrap(qnorm(as<NumericVector>(wrap(dat)))));
    xev.reshape(m, d);
    xdat.reshape(n, d);
    
    // apply bandwidth matrix
    mat zev  = (inv(B) * (xev).t()).t();
    mat zdat = (inv(B) * (xdat).t()).t();
    
    // create temporaray objects for loop
    mat tmpmat(n, d);
    rowvec xevi;
    double tmp;
    NumericVector rescale(2);
    
    for(int i = 0; i < m; ++i){
        // compute standard kernel estimate
        tmpmat = zdat - repmat(zev.row(i), n, 1);
        tmp = mean(kern_gauss(tmpmat, rep(1.0, 2)));
        
        // rescale to obtain a kernel estimate of copula density
        xevi = xev.row(i);
        rescale = dnorm(as<NumericVector>(wrap(xevi)));
        out[i] = tmp / (rescale[0] * rescale[1] * det(B));
    }
    
    return out;
}



