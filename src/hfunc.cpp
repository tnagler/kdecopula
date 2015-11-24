#include <RcppArmadillo.h>
#include <interp.h>
#include <integrate.h>
#include <hfunc.h>
using namespace Rcpp;

//// evaluate h-function via integration
// [[Rcpp::export]]
Rcpp::NumericVector eval_hfunc_2d(const Rcpp::NumericMatrix& uev,
                                  const int& cond_var, 
                                  const Rcpp::NumericMatrix& vals,
                                  const Rcpp::NumericVector& grid) 
{
    int N = uev.nrow();
    int m = grid.size();
    NumericVector tmpvals(m), out(N);
    NumericMatrix tmpgrid(m, 2);
    double upr = 0.0;
    double tmpint, int1; 
    tmpint = 0.0;
    
    for (int n = 0; n < N; ++n) {
        if (cond_var == 1) {
            upr = uev(n, 1);
            tmpgrid(_, 0) = rep(uev(n, 0), m);
            tmpgrid(_, 1) = grid;
        } else if (cond_var == 2) {
            upr = uev(n, 0);
            tmpgrid(_, 0) = grid;
            tmpgrid(_, 1) = rep(uev(n, 1), m);
        }
        tmpvals = interp_2d(tmpgrid, vals, grid);
        tmpint = int_on_grid(upr, tmpvals, grid);
        int1 = int_on_grid(1.0, tmpvals, grid);
        out[n] = tmpint/int1;
        out[n] = fmax(out[n], 1e-10);
        out[n] = fmin(out[n], 1-1e-10);
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector eval_hfunc(const Rcpp::NumericMatrix& uev, 
                               const Rcpp::IntegerVector& cond_var,
                               const Rcpp::IntegerVector& uncond_var,
                               const Rcpp::NumericVector& vals,
                               const Rcpp::NumericVector& grid, 
                               const Rcpp::NumericMatrix& helpgrid,
                               const Rcpp::IntegerMatrix& helpind)
{
    int N = uev.nrow();
    int m = grid.size();
    IntegerVector gridsq = seq_len(m) - 1;
    NumericVector newvals(helpgrid.nrow()), newvals1(helpgrid.nrow());
    NumericVector tmpvals(m), tmpvals1(m);
    NumericMatrix uev1 = clone(uev);
    NumericVector out(N);
    NumericMatrix hg_new(clone(helpgrid));
    int uncond_len = uncond_var.size();
    
    for (int n = 0; n < N; ++n) {
        // insert the conditioning values in helpgrid 
        for (int i = 0; i < cond_var.size(); ++i) {
            hg_new(_, cond_var[i]-1) = rep(uev(n, cond_var[i]-1), helpgrid.nrow());
        }
        
        // interpolate on helpgrid
        newvals = interp(hg_new, vals, grid, helpind);
        newvals1 = clone(newvals);
        
        // recursively integrate over a sequence of m values
        for (int j = 0; j < uncond_len; ++j){
            for (int p = 0; p < pow((double)m, uncond_len - j - 1); ++p) {
                tmpvals = newvals[p*m + gridsq];
                tmpvals1 = newvals1[p*m + gridsq];
                newvals[p] = int_on_grid(uev(n, uncond_var[j]-1), tmpvals, grid);
                newvals1[p] = int_on_grid(1.0, tmpvals1, grid);
            }
        }
        out[n] = newvals[0]/newvals1[0];
        out[n] = fmax(out[n], 1e-10);
        out[n] = fmin(out[n], 1-1e-10);
    }
    return out;
}


//// invert h-function 
// [[Rcpp::export]]
Rcpp::NumericVector inv_hfunc(const Rcpp::NumericMatrix& uev,
                              const int& cond_var,
                              const Rcpp::NumericMatrix& vals,
                              const Rcpp::NumericVector& grid) 
{
    int N = uev.nrow();
    int m = grid.size();
    NumericVector tmpvals(m), out(N);
    NumericMatrix tmpgrid(m, 2);
    double upr = 0.0;
    
    for (int i = 0; i < N; ++i) {
        if (cond_var == 1) {
            upr = uev(i, 1);
            tmpgrid(_, 0) = rep(uev(i, 0), m);
            tmpgrid(_, 1) = grid;
            
        } else if (cond_var == 2) {
            upr = uev(i, 0);
            tmpgrid(_, 0) = grid;
            tmpgrid(_, 1) = rep(uev(i, 1), m);
        }
        tmpvals = interp_2d(tmpgrid, vals, grid);
        out[i] = inv_int_on_grid(upr, tmpvals, grid);
        out[i] = fmax(out[i], 1e-10);
        out[i] = fmin(out[i], 1-1e-10);
    }
    return out;
}


// full d-dimensional integrals (copula cdf)
// [[Rcpp::export]]
Rcpp::NumericVector eval_cdf(const Rcpp::NumericMatrix& uev, 
                             const Rcpp::NumericVector& vals,
                             const Rcpp::NumericVector& grid, 
                             const Rcpp::NumericMatrix& helpgrid, 
                             const Rcpp::IntegerMatrix& helpind)
{
    int N = uev.nrow();
    int d = uev.ncol();
    int m = grid.size();
    IntegerVector gridsq = seq_len(m) - 1;
    
    NumericVector newvals(helpgrid.nrow()), newvals1(helpgrid.nrow());
    NumericVector tmpvals(m), tmpvals1(m);
    NumericVector out(N);
    
    for (int n = 0; n < N; ++n) {
        // interpolate on helpgrid
        newvals = interp(helpgrid, vals, grid, helpind);
        newvals1 = clone(newvals); /* also compute integral to one for rescaling */
        // recursively integrate over a sequence of m values
        for (int j = 0; j < d; ++j){
            for (int p = 0; p < pow((double)m, d - j - 1); ++p) {
                tmpvals = newvals[p*m + gridsq];
                tmpvals1 = newvals1[p*m + gridsq];
                newvals[p] = int_on_grid(uev(n, j), tmpvals, grid);
                newvals1[p] = int_on_grid(1.0, tmpvals1, grid);
            }
        }
        
        out[n] = newvals[0]/newvals1[0];
        out[n] = fmax(out[n], 1e-10);
        out[n] = fmin(out[n], 1-1e-10);
    }
    return out;
}
