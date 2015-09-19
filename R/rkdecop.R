#' Simulate data from a kdecopula object
#' 
#' @param n integer; number of observations.
#' @param quasi logical; the default (\code{FALSE}) returns pseudo-random
#' numbers, use \code{TRUE} for quasi-random numbers (generalized Halton, see
#' \code{\link[qrng:ghalton]{ghalton}}).
#' 
#' @rdname kdecop
#' 
rkdecop <- function(n, obj, quasi = FALSE) {
    n <- round(n)
    stopifnot(inherits(obj, "kdecopula"))
    stopifnot(is.logical(quasi))
    
    if (!quasi) {
        # simulate independent uniform random variables
        W <- cbind(runif(n), runif(n))
    } else {
        # generate quasi random numbers
        w <- ghalton(n, d = 2)
    }
    
    # if independence copula is specified, return W
    if("indep.copula" %in% class(obj))
        return(W)
    
    # invert h-function otherwise
    U2 <- inv_hfunc(W,
                    1L,
                    obj$estimate,
                    obj$grid)
    
    ## return results
    out <- cbind(W[, 1], U2)
    colnames(out) <- NULL
    out
}
