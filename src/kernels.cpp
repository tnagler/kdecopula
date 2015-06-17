#include <RcppArmadillo.h>
#include <kernels.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector kern_epan(NumericVector x){
    NumericVector out(x.size());
    for(int i = 0; i < x.size(); ++i){
        if((std::fabs(x[i]) >= 1.0)){
            out[i] = 0;
        } else{
            out[i] = 3.0/4.0 * (1 - pow(x[i], 2)); 
        } 
    }
    return out;
}

// [[Rcpp::export]]
NumericVector kern_epan_2d(NumericVector x, NumericVector y){
    return(kern_epan(x) * kern_epan(y));
}

// [[Rcpp::export]]
NumericVector kern_gauss(NumericVector x){
    NumericVector out(x.size());
    for(int i = 0; i < x.size(); ++i){
        if((std::fabs(x[i]) >= 5)){
            out[i] = 0;
        } else{
            out[i] = exp(- 0.5 * pow(x[i], 2)) / (sqrt(2 * M_PI)) / 0.9999994267; 
        } 
    }
    return out;
}

// [[Rcpp::export]]
NumericVector kern_gauss_2d(NumericVector x, NumericVector y){
    return(kern_gauss(x) * kern_gauss(y));
}
