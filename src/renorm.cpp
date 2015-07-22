#include <RcppArmadillo.h>
#include <interp.h>
#include <integrate.h>
#include <renorm.h>
using namespace Rcpp;

//// renormalize a subset of an array
// [[Rcpp::export]]
NumericVector ren_subs(NumericVector vals, NumericVector grid, int d) {
    int m = grid.size();
    int N = vals.size();
    IntegerVector gridsq = seq_len(m) - 1;
    NumericVector newvals = clone(vals);
    NumericVector tmpvals(m);
    NumericVector out(N);
    
    // recursively integrate over a sequence of m values
    for (int j = 0; j < d; ++j){
        for (int p = 0; p < pow((double)m, d - j - 1); ++p) {
            tmpvals = newvals[p*m + gridsq];
            newvals[p] = int_on_grid(1.0, tmpvals, grid);
        }
    }
    
    for (int i = 0; i < N; ++i) 
    out[i] = vals[i] / fmax(newvals[0], 1e-10);
    
    return out;
}

//// rescale d-dimensional copula density
// the matrix 'helpind' is the (d-1)-dimensional expanded sequence 0:(knots-1) 
// [[Rcpp::export]]
NumericVector renorm(NumericVector x, NumericVector grid, int times, IntegerMatrix helpind) {
    double m = grid.size();
    IntegerVector dims = x.attr("dim");
    int d = dims.size();
    int mtd = pow(m, d-1);
    
    NumericVector M(x.size()), vals;
    M = clone(x);
    IntegerMatrix tmpind(mtd, d);
    
    for(int t = 0; t < times; ++t) {
        for (int j = 0; j < d; ++j) {
            for (int k = 0; k < m; ++k) {
                // get indices for subsetting array M
                int pl = 0;
                for (int jj = 0; jj < d; ++jj) {
                    if (jj == j) {
                        tmpind(_, jj) = rep(k, mtd);
                    } else { 
                        tmpind(_, jj) = helpind(_, pl);
                        ++pl; 
                    }
                }
                vals = M[get(tmpind, dims)];
                M[get(tmpind, dims)] = ren_subs(vals, grid, d-1);
            }
        }
    }
    return M;
}

//// calculate errors for uniform margins
NumericMatrix calc_errors(NumericVector x, NumericVector grid, int times, IntegerMatrix helpind) {
    double m = grid.size();
    IntegerVector dims = x.attr("dim");
    int d = dims.size();
    int mtd = pow(m, d-1);
    NumericMatrix err(times, m);
    
    NumericVector M(x.size());
    M = clone(x);
    IntegerMatrix tmpind(mtd, d);
    
    for(int t = 0; t < times; ++t) {
        for (int j = 0; j < d; ++j) {
            for (int k = 0; k < m; ++k) {
                // get indices for subsetting array M
                int pl = 0;
                for (int jj = 0; jj < d; ++jj) {
                    if (jj == j) {
                        tmpind(_, jj) = rep(k, mtd);
                    } else { 
                        tmpind(_, jj) = helpind(_, pl);
                        ++pl; 
                    }
                }
                if (j == 0) {
                    err(t, k) = int_on_grid(1.0, M[get(tmpind, dims)], grid);
                }
                M[get(tmpind, dims)] = ren_subs(M[get(tmpind, dims)], grid, d-1);
            }
        }
    }
    return err;
}

