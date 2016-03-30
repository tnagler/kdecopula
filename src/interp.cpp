#include <RcppArmadillo.h>
#include <interp.h>
using namespace Rcpp;


//// evaluate a cubic polynomial with given coefficents
double cubic_poly(const double& x, 
                  const Rcpp::NumericVector& a) 
{
    double x2 = x*x;
    double x3 = x2*x;
    return a[0] + a[1]*x + a[2]*x2 + a[3]*x3;
}

//// indefinite integral of cubic polynomial
double cubic_indef_integral(const double& x,
                            const Rcpp::NumericVector& a) 
{
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    return a[0]*x + a[1]/2.0*x2 + a[2]/3.0*x3 + a[3]/4.0*x4;
}

//// integral of cubic polynomial from lowr to upr
double cubic_integral(const double& lowr,
                      const double& upr, 
                      const Rcpp::NumericVector& a) 
{
    return cubic_indef_integral(upr, a) - cubic_indef_integral(lowr, a);
}

//// numerically invert a cubic integral (with 0 as lower bound)
// -> maximum of 500 search iterations
double inv_cubic_integral(const double& q,
                          const Rcpp::NumericVector& a) 
// {
//     double x = fabs(q);
//     double qtest, tmpint, sign;
//     for (int i = 0; i < 100; ++i) {
//         tmpint = cubic_integral(0.0, x, a);
//         if (tmpint > 0) sign = 1; else sign = -1;
//         qtest = fabs(tmpint) - q;
//         if (fabs(qtest) < 1e-20) break;
//         x += - sign * qtest / cubic_poly(x, a);
//     }
//     return x;
// }
{
    double x0, x1, ql, qh, ans, val;
    ans = 0.0, val = 0.0; x0 = 0.0; x1 = 1.0;
    ql = 0.0;
    qh = cubic_integral(0.0, x1, a);
    int br = 0;
    double tol = ::fmax(1e-10 * (x1 - x0), 1e-10);
    
    // check if already close enough (or 1.0 is exceeded)
    ql = ql - q;
    qh = qh - q;
    if (::fabs(ql) <= tol) {
        ans = x0;
        br = 1;
    }
    if ((::fabs(qh) <= tol) | (qh < 0)) {
        ans = x1;
        br = 1;
    }
    
    for (int it = 0; it < 20; ++it) {
        ans = (x0 + x1) / 2.0;
        val = cubic_integral(0.0, ans, a);
        val = val - q;
        // stop if values become too close (avoid infinite loop)
        if (::fabs(val) <= tol)
            br = 1;
        if (::fabs(x0 - x1) <= tol)
            br = 1;
        // check which of x0, x1 is closer
        if (val > 0.0) {
            x1 = ans;
            qh = val;
        } else {
            x0 = ans;
            ql = val;
        }
        // stop if convergence 
        if (br == 1)
            break;
    }
    
    return ans;
}


//// calculate coefficients for cubic splines (input must have length 4)
Rcpp::NumericVector coef(const Rcpp::NumericVector& vals,
                         const Rcpp::NumericVector& grid) 
{
    NumericVector a(4);
    
    double dt0 = grid[1] - grid[0];
    double dt1 = grid[2] - grid[1];
    double dt2 = grid[3] - grid[2];
    
    /* check for repeated points (important for boundaries) */
    if (dt1 < 1e-4) dt1 = 1.0;
    if (dt0 < 1e-4) dt0 = dt1;
    if (dt2 < 1e-4) dt2 = dt1;
    
    // compute tangents when parameterized in [t1,t2]
    double dx1 = (vals[1] - vals[0]) / dt0 - (vals[2] - vals[0]) / (dt0 + dt1) + (vals[2] - vals[1]) / dt1;
    double dx2 = (vals[2] - vals[1]) / dt1 - (vals[3] - vals[1]) / (dt1 + dt2) + (vals[3] - vals[2]) / dt2;
    
    // rescale tangents for parametrization in [0,1]
    dx1 *= dt1;
    dx2 *= dt1;
    
    // compute coefficents
    a[0] = vals[1];
    a[1] = dx1;
    a[2] = -3*vals[1] + 3*vals[2] - 2*dx1 - dx2;
    a[3] = 2*vals[1] - 2*vals[2] + dx1 + dx2;
    
    return a;
} 

//// interpolate in one dimension (inputs must have length 4)
double interp_on_grid(const double& x,
                      const Rcpp::NumericVector& vals, 
                      const Rcpp::NumericVector& grid) {
    NumericVector a = coef(vals, grid);
    double xev = fmax((x - grid[1]), 0)/(grid[2] - grid[1]);
    return cubic_poly(xev, a);
}

//// interpolate in two dimensions
// [[Rcpp::export]]
Rcpp::NumericVector interp_2d(const Rcpp::NumericMatrix& x, 
                              const Rcpp::NumericMatrix& vals, 
                              const Rcpp::NumericVector& grid)
{
    int N = x.nrow();
    int m = grid.size();
    NumericVector y(4), tmpvals(4), tmpgrid(4), out(N);
    int i = 0;
    int j = 0;
    int i0, i2;
    
    for (int n = 0; n < N; ++n) {
        // find cell
        for (int k = 1; k < (m-1); ++k) {
            if((x(n, 0) >= grid[k])) i = k;
            if((x(n, 1) >= grid[k])) j = k;
        }    
        
        // construct grid for first direction
        i0 = std::max(i-1, 0);
        i2 = std::min(i+2, m-1);
        tmpgrid = NumericVector::create(grid[i0], grid[i], grid[i+1], grid[i2]);
        
        // interpolate in one direction (four times)
        for(int s = 0; s < 4; ++s) {
            i0 = std::max(i-1, 0);
            i2 = std::min(i+2, m-1);
            int jj = std::min(m-1, j-1+s);
            jj = std::max(0, jj);
            
            tmpvals = NumericVector::create(vals(i0, jj),
                                            vals(i,   jj),
                                            vals(i+1, jj),
                                            vals(i2,  jj));
            
            y[s] = interp_on_grid(x(n, 0), tmpvals, tmpgrid);
            y[s] = fmax(y[s], 0.0);
        }
        
        // use these four points to interpolate in the remaining direction#
        i0 = std::max(j-1, 0);
        i2 = std::min(j+2, m-1);
        tmpgrid = NumericVector::create(grid[i0], grid[j], grid[j+1], grid[i2]);
        out[n] = interp_on_grid(x(n, 1), y, tmpgrid);
        out[n] = fmax(out[n], 1e-15);
    }
    
    return out;
}

// subset array in C++
Rcpp::IntegerVector get(const Rcpp::IntegerMatrix& ind, 
                        const Rcpp::IntegerVector& dims) 
{
    int ndims = dims.size();
    int m = dims[0];
    int N = ind.nrow();
    IntegerVector tmpi(ndims);
    IntegerVector out(N);
    
    // calcualte correct index for vector instead of array
    for (int n = 0; n < N; ++n) {
        for (int i = 0; i < ndims; ++i) {
            tmpi[i] = ind(n, i) * pow((double)m, i);
        }
        out[n] = sum(tmpi);
    }
    
    return out;
} 

//// interpolate on an array with more than two dimensions
// [[Rcpp::export]]
Rcpp::NumericVector interp(const Rcpp::NumericMatrix& x,
                           const Rcpp::NumericVector& vals, 
                           const Rcpp::NumericVector& grid, 
                           const Rcpp::IntegerMatrix& helpind) 
{
    int d = x.ncol();
    int N = x.nrow();
    int m = grid.size();
    IntegerVector dims = vals.attr("dim");
    NumericVector y(4), tmpvals(4), tmpgrid(4), out(N);
    IntegerMatrix tmpi(4, d); 
    int i0, i2;
    int ii;
    
    IntegerVector i = rep(0, d);
    for (int n = 0; n < N; ++n) { 
        // find cell
        for (int k = 2; k < (m-1); ++k) {
            for (int j = 0; j < d; ++j) {
                if ((x(n, j) >= grid[k]))
                    i[j] = k;
            }
        } 
        
        // store grid for first  direction
        NumericVector newvals(helpind.nrow()/4);
        i0 = std::max(i[0]-1, 0);
        i2 = std::min(i[0]+2, m-1);
        tmpgrid = NumericVector::create(grid[i0], grid[i[0]], grid[i[0]+1], grid[i2]); 
        // get values and interpolate
        for (int p = 0; p < helpind.nrow()/4; ++p) {
            IntegerMatrix tmpi(4, d);
            for (int kk = 0; kk < 4; ++kk) {
                for (int jj = 0; jj < d; ++jj){
                    ii = i[jj] + helpind(p*4 + kk, jj);
                    ii = std::max(ii, 0);
                    ii = std::min(ii, m-1);
                    tmpi(kk, jj) = ii;
                }        
            }   
            tmpvals = vals[get(tmpi, dims)];
            newvals[p] = interp_on_grid(x(n, 0), tmpvals, tmpgrid);
            newvals[p] = fmax(newvals[p], 0.0);
        }
        
        // do the same for remain directions
        for (int j = 1; j < d; ++j){
            i0 = std::max(i[j]-1, 0);
            i2 = std::min(i[j]+2, m-1);
            tmpgrid = NumericVector::create(grid[i0], grid[i[j]], grid[i[j]+1], grid[i2]);
            for (int p = 0; p < helpind.nrow()/pow(4.0, j + 1); ++p) {
                tmpvals = NumericVector::create(newvals[p*4],
                                                newvals[p*4 + 1],
                                                       newvals[p*4 + 2],
                                                              newvals[p*4 + 3]);
                newvals[p] = interp_on_grid(x(n, j), tmpvals, tmpgrid);
                newvals[p] = fmax(newvals[p], 0.0);
            }
        }
        
        
        out[n] = newvals[0];
        out[n] = fmax(out[n], 1e-15);
    } 
    
    return out;
}

