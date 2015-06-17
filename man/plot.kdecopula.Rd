\name{plot.kdecopula}
\alias{plot.kdecopula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting \code{kdecopula} objects
}
\description{
Produces contour or perspective plots for a \code{kdecopula} object.
}
\usage{
\method{plot}{kdecopula}(x, type = "contour", margins, size, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{kdecopula} object.
}
\item{type}{
plot type; either \code{"contour"} or \code{"surface"} (partial matching is activated) for a contour or perspective/surface plot respectively.
}
  \item{margins}{ either \code{"unif"} for the original copula density or \code{"norm"} for the transformed density with standard normal margins (partial matching is activated). Default is \code{"norm"} for \code{type = "contour"}, and \code{"unif"} for \code{type = "surface"}.
} 
  \item{size}{
integer; the plot is based on values on a \eqn{size x size} grid; default is 100 for \code{type = "contour"}, and 25 for \code{type = "surface"}.
}
  \item{\dots}{
optional arguments passed to \code{\link{contour}} or \code{\link{wireframe}}.}
}
\author{
Thomas Nagler
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[kdecopula:kdecop]{kdecop}},
\code{\link[graphics:contour]{contour}},
\code{\link[lattice:wireframe]{wireframe}}

}
\examples{
## load data and transform with empirical cdf
data(wdbc)
udat <- apply(wdbc[, -1], 2, function(x) rank(x)/(length(x)+1))

## estimation of copula density of variables 5 and 6
obj <- kdecop(udat[, 5:6])

## plots
plot(obj)  # (marginal normal) contour plot
plot(obj, margins = "unif")  # contour plot of actual copula density
plot(obj, type = "surf")  # surface plot of actual copula density 


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{Copula density estimation}
