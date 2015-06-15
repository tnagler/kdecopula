#include <RcppArmadillo.h>
#include <kernels.h>
#include <interp.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
NumericVector eval_mr(NumericMatrix uev, NumericMatrix dat, double b) {
    double n = dat.nrow();
    NumericVector out(uev.nrow());
    NumericVector tmp1(n), tmp2(n), tmp3(n),
    tmp4(n), tmp5(n), tmp6(n),
    tmp7(n), tmp8(n), tmp9(n);
    
    for (int i = 0; i < uev.nrow(); ++i) {
        tmp1 = kern_epan_2d((uev(i, 0) - dat(_, 0)) / b, (uev(i, 1) - dat(_, 1)) / b); 
        
        tmp2 = kern_epan_2d((uev(i, 0) + dat(_, 0)) / b, (uev(i, 1) - dat(_, 1)) / b);
        tmp3 = kern_epan_2d((uev(i, 0) - dat(_, 0)) / b, (uev(i, 1) + dat(_, 1)) / b);
        tmp4 = kern_epan_2d((uev(i, 0) + dat(_, 0)) / b, (uev(i, 1) + dat(_, 1)) / b);
        
        tmp5 = kern_epan_2d((uev(i, 0) - dat(_, 0)) / b, (uev(i, 1) + dat(_, 1) - 2) / b);
        tmp6 = kern_epan_2d((uev(i, 0) + dat(_, 0)) / b, (uev(i, 1) + dat(_, 1) - 2) / b);
        
        tmp7 = kern_epan_2d((uev(i, 0) + dat(_, 0) - 2) / b, (uev(i, 1) - dat(_, 1)) / b);
        tmp8 = kern_epan_2d((uev(i, 0) + dat(_, 0) - 2) / b, (uev(i, 1) + dat(_, 1)) / b);
        
        tmp9 = kern_epan_2d((uev(i, 0) + dat(_, 0) - 2) / b, (uev(i, 1) + dat(_, 1) - 2) / b);
        
        out[i] = sum(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8) / (pow(b, 2) * n);
    }
    return out;
}

// [[Rcpp::export]]
NumericVector eval_beta(NumericMatrix uev, NumericMatrix dat, double b){
    NumericVector out(uev.nrow());
    double n = dat.nrow();
    NumericVector B0(n), B1(n);
    
    for (int i = 0; i < uev.nrow(); ++i) {
        NumericVector B = rep(1.0, n);
        for (int j = 0; j < uev.ncol(); ++j) {
            B = B * dbeta(dat(_, j), uev(i, j) / b + 1, (1 - uev(i, j)) / b + 1);
        }
        out[i] = sum(B) / n;
    }
    return out;
}


// [[Rcpp::export]]
NumericVector eval_trafo(NumericMatrix uev, NumericMatrix dat, arma::mat B){
    NumericVector out(uev.nrow());
    double n = dat.nrow();
    int d = uev.ncol();
    NumericMatrix xev(uev.nrow(), uev.ncol());
    NumericMatrix xdat(dat.nrow(), dat.ncol());
    
    for(int j = 0; j < d; ++j){
        xev(_, j) = qnorm(as<NumericVector>(wrap(uev(_, j))));
        xdat(_, j) = qnorm(as<NumericVector>(wrap(dat(_, j))));
    }
    
    
    NumericMatrix zev  = wrap((inv(B) * as<arma::mat>(xev).t()).t());
    NumericMatrix zdat = wrap((inv(B) * as<arma::mat>(xdat).t()).t());
    
    NumericMatrix tmpmat(dat.nrow(), dat.ncol());
    double tmp;
    NumericVector rescale(2);
    for(int i = 0; i < uev.nrow(); ++i){
        for(int j = 0; j < 2; ++j){
            tmpmat(_, j) = zev(i, j) - zdat(_, j);
        }
        tmp = sum(kern_gauss_2d(tmpmat(_, 0), tmpmat(_, 1))) / n;
        rescale = dnorm(as<NumericVector>(wrap(xev(i, _))));
        out[i] = tmp / (rescale[0] * rescale[1] * det(B));
    }
    
    return out;
}



