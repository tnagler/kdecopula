\name{kdecopula-package}
\alias{kdecopula-package}
\docType{package}
\title{
Kernel Smoothing for Bivariate Copula Densities
}
\description{
The package provides fast implementations of kernel smoothing techniques for bivariate copula densities, in particular density estimation and resampling.
}
\details{
\tabular{ll}{
Package: \tab kdecopula\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2015-06-15\cr
License: \tab GPL-3\cr
Maintainer: \tab Thomas Nagler <thomas.nagler@tum.de>\cr
Depends: \tab R (>= 3.0.0)\cr
Imports: \tab lattice, cubature, locfit, Rcpp (>= 0.11.2), VineCopula (>= 1.2)\cr
LinkingTo: \tab Rcpp, RcppArmadillo\cr
}
The function \code{\link{kdecop}} can be used to estimate a copula density with a number of popular kernel estimators. The density estimate can be evaluated on arbitrary points with \code{\link[kdecopula:dkdecop]{dkdecop}}; the cdf with \code{\link[kdecopula:pkdecop]{pkdecop}}. Furthermore, synthetic data can be simulated with \code{\link[kdecopula:rkdecop]{rkdecop}}, and several plot options are provided by
\code{\link[kdecopula:plot.kdecopula]{plot.kdecopula}}.
}
\author{
Thomas Nagler
}
\references{
Gijbels, I. and Mielniczuk, J. (1990). Estimating the density of a copula function.
Communications in Statistics - Theory and Methods, 19(2):445-464. \cr 
\cr
Charpentier, A., Fermanian, J.-D., and Scaillet, O. (2006). The estimation of copulas:
Theory and practice. In Rank, J., editor, Copulas: From theory to application in
finance. Risk Books. \cr 
\cr
Geenens, G., Charpentier, A., and Paindaveine, D. (2014). Probit transformation for
nonparametric kernel estimation of the copula density. arXiv:1404.4414 [stat.ME]. \cr
\cr
Nagler, T. (2014). "Kernel Methods for Vine Copula Estimation", Master's Thesis, Technische Universitaet Muenchen, \url{https://mediatum.ub.tum.de/node?id=1231221}
}
\keyword{ package }
\examples{
## load data and transform with empirical cdf
data(wdbc)
udat <- apply(wdbc[, -1], 2, function(x) rank(x)/(length(x)+1))

## estimation of copula density of variables 5 and 6
dens.est <- kdecop(udat[, 5:6])
plot(dens.est) 

## evaluate density estimate at (u1,u2)=(0.123,0.321)
dkdecop(c(0.123, 0.321), dens.est) 

## evaluate cdf estimate at (u1,u2)=(0.123,0.321)
pkdecop(c(0.123, 0.321), dens.est) 

## simulate 500 samples from density estimate
rkdecop(500, dens.est)
}
