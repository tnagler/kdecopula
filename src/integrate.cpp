#include <RcppArmadillo.h>
#include <interp.h>
#include <integrate.h>
using namespace Rcpp;


//' Integrate a spline interpolant
//' 
//' @param upr upper limit of integration (lower is 0).
//' @param vals vector of values to be interpolated and integrated.
//' @param grid vector of grid points on which vals has been computed.
//' 
//' @return Integral of interpolation spline defined by (vals, grid).
//' 
//' @noRd
// [[Rcpp::export]]
double int_on_grid(const double& upr, 
                   const Rcpp::NumericVector& vals,
                   const Rcpp::NumericVector& grid) 
{
    int m = grid.size();
    NumericVector tmpvals(4), tmpgrid(4), tmpa(4), a(4);
    double uprnew, newint;
    
    double tmpint = 0.0;
    
    // linear extrapolation for interval [0, grid[0]]
    double dx = (vals[1] - vals[0]) / (grid[1] - grid[0]);
    dx *= (grid[1] - grid[0]);
    double firstlen = fmin(grid[0], upr);
    tmpint += vals[0]*firstlen + (dx/2.0 * pow(firstlen, 2.0))*grid[0];
    
    if (upr > grid[0]) {
        // go up the grid and integrate
        for (int k = 0; k < m-1; ++k) {
            // stop loop if fully integrated
            if (upr < grid[k]) break;
            
            // select length 4 subvectors and calculate spline coefficients
            tmpvals[0] = vals[std::max(k-1, 0)]; 
            tmpvals[1] = vals[k];
            tmpvals[2] = vals[k+1];
            tmpvals[3] = vals[std::min(k+2, m-1)];
            
            tmpgrid[0] = grid[std::max(k-1, 0)]; 
            tmpgrid[1] = grid[k];
            tmpgrid[2] = grid[k+1];
            tmpgrid[3] = grid[std::min(k+2, m-1)];
            
            tmpa = coef(tmpvals, tmpgrid, a);
            
            // don't integrate over full cell if upr is in interior
            uprnew = (upr - grid[k]) / (grid[k+1] - grid[k]);
            newint = cubic_integral(0.0, fmin(1.0, uprnew), tmpa);
            tmpint += newint * (grid[k+1] - grid[k]);
        }
        
        // if grid is not sufficient, extrapolate linearly
        if(upr > grid[m-1]) {        
            double dx = (vals[m-1] - vals[m-2]) / (grid[m-1] - grid[m-2]);
            dx *= (grid[m-1] - grid[m-2]);  
            uprnew = (upr - grid[m-1]) / (1.0 - grid[m-2]);
            double lastlen = upr - grid[m-1];
            tmpint += vals[m-1]*lastlen + (dx/2.0 * pow(uprnew, 2.0))*(1.0 - grid[m-1]);
        } 
    } 
    
    return tmpint;
}

// The following was used in previous versions but is deprecated for now.
//
// //' Inverse of the integral of a spline interpolant
// //' 
// //' @param qq argument of the inverse integral (the 'quantile').
// //' @param vals vector of values to be interpolated and integrated.
// //' @param grid vector of grid points on which vals has been computed.
// //' 
// //' @return Integral of interpolation spline defined by (vals, grid).
// //' 
// //' @noRd
// // [[Rcpp::export]]
// double inv_int_on_grid(const double& qq, 
//                        const Rcpp::NumericVector& vals,
//                        const Rcpp::NumericVector& grid)
// {
//     int m = grid.size();
//     NumericVector tmpvals(4), tmpgrid(4), tmpa(4), a(4);
//     double uprnew, newint, out, qtest;
//     double tmpint = 0.0;
//     int tmpk = 0;
//     double q = qq;
//     
//     q *= int_on_grid(1.0, vals, grid);
//     
//     double dx = (vals[1] - vals[0]) / (grid[1] - grid[0]);
//     dx *= (grid[1] - grid[0]);
//     tmpint += vals[0]*grid[0] + (dx/2.0 * pow(grid[0], 2.0))*grid[0];
//     qtest = tmpint;
//     uprnew = (q - qtest)/(grid[2] - grid[1]);
//     
//     // go up the grid and integrate as long as target value is above integral value
//     if (q > qtest) {
//         for (int k = 1; k < m-2; ++k) {
//             // select length 4 subvectors and calculate spline coefficients
//             tmpvals[0] = vals[k-1]; 
//             tmpvals[1] = vals[k];
//             tmpvals[2] = vals[k+1];
//             tmpvals[3] = vals[k+2];
//             
//             tmpgrid[0] = grid[k-1]; 
//             tmpgrid[1] = grid[k];
//             tmpgrid[2] = grid[k+1];
//             tmpgrid[3] = grid[k+2];
//             
//             tmpa = coef(tmpvals, tmpgrid, a);
//             newint = cubic_integral(0.0, 1.0, tmpa);
//             tmpint += newint * (grid[k+1] - grid[k]);
//             tmpk = k;
//             if (tmpint > q) break;
//             uprnew = (q - tmpint)/(grid[tmpk+1] - grid[tmpk]);
//         }
//     } else {
//         return 2.0/dx*sqrt(q) * grid[0];
//     }
//     
//     // solve in cell
//     double lastgrid;
//     if (tmpint > q) {
//         lastgrid = grid[tmpk];
//         out = lastgrid + inv_cubic_integral(uprnew, tmpa)*(grid[tmpk+1] - grid[tmpk]);    
//     } else {
//         double dx = (vals[m-1] - vals[m-2]) / (grid[m-1] - grid[m-2]);
//         dx *= (grid[m-1] - grid[m-2]);  
//         out = grid[tmpk+1] + 2.0/dx*sqrt(uprnew)*(1.0 - grid[m-1]);
//     }
//     return out;
// }
