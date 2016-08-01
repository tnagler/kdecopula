#include <RcppArmadillo.h>
#include <kernels.h>
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
    NumericVector tmpvals(m), out(N), tmpa(4), tmpb(4);
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
        tmpvals = interp_2d(tmpgrid, vals, grid, tmpa, tmpb);
        tmpint = int_on_grid(upr, tmpvals, grid);
        int1 =  int_on_grid(1.0, tmpvals, grid);
        out[n] = tmpint/int1;
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
    NumericVector out(uev.nrow()), ans, val, x0, x1;
    double q;
    int m = grid.size();
    NumericMatrix tmpu0(1, 2), tmpu1(1, 2);
    ans = 0.0, val = 0.0;
    double tol = 1e-10;
    
    for (int i = 0; i < uev.nrow(); ++i) {
        int br = 0;
        x0 = 0;
        x1 = 1;
        
        // evaluation points at boundary
        if (cond_var == 1) {
            q = uev(i, 1);
            tmpu0(0, 0) = uev(i, 0);
            tmpu0(0, 1) = x0[0];
            tmpu1(0, 0) = uev(i, 0);
            tmpu1(0, 1) = x1[0];
            
        } else if (cond_var == 2) {
            q = uev(i, 0);
            tmpu0(0, 0) = x0[0];
            tmpu0(0, 1) = uev(i, 1);
            tmpu1(0, 0) = x1[0];
            tmpu1(0, 1) = uev(i, 1);
        }
        
        // evaluate h-function at boundary points
        NumericVector ql = eval_hfunc_2d(tmpu0, cond_var, vals, grid);
        NumericVector qh = eval_hfunc_2d(tmpu1, cond_var, vals, grid);
        ql = ql - q;
        qh = qh - q;
        
        // check if already close enough (unless at boundary)
        if ((::fabs(ql[0]) < tol) && (q > 1e-9)) {
            ans = x0;
            br = 1;
        } else if ((::fabs(qh[0]) < tol) && (q < 1-1e-9)) {
            ans = x1;
            br = 1;
        }
        
        // find inverse by bisection
        int maxit = 15;
        for (int it = 0; it < maxit; ++it) {
            // set new evaluation point
            ans[0] = (x0[0] + x1[0]) / 2.0;
            if (cond_var == 1) {
                q = uev(i, 1);
                tmpu0(0, 0) = uev(i, 0);
                tmpu0(0, 1) = ans[0];
                
            } else if (cond_var == 2) {
                q = uev(i, 0);
                tmpu0(0, 0) = ans[0];
                tmpu0(0, 1) = uev(i, 1);
            }
            
            // evaluate h-function
            NumericVector val = eval_hfunc_2d(tmpu0, cond_var, vals, grid);
            val[0] = val[0] - q;

            // find section for next iteration
            if (::fabs(val[0]) < 1e-9) {
                if (q <= 9e-9) {
                    // go to upper section if q == 1e-10 
                    x0[0] = ans[0];
                    ql[0] = val[0];
                } else if (q >= 1 - 9e-9)  {
                    // go to lower section if q == 1 - 1e-10 
                    x1[0] = ans[0];
                    qh[0] = val[0];
                } else {
                    br = 1;
                }
            } else if (val[0] > 0.0) {
                x1[0] = ans[0];
                qh[0] = val[0];
            } else if (val[0] < 0.0) {
                x0[0] = ans[0];
                ql[0] = val[0];
            } 
            
            // stop if values become too close
            if (::fabs(x0[0] - x1[0]) <= tol)
                br = 1;
            
            if (br == 1)
                break;
        }
        
        out[i] = ans[0];
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
