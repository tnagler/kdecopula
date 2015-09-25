#' Evaluate the density of a kdecopula object
#' 
#' @param stable logical; option for stabilizing the estimator: the estimated
#' density is cut off at \eqn{50}.
#' @rdname kdecop
#' 
dkdecop <- function(u, obj, stable = FALSE) {
    stopifnot(is.numeric(u))
    stopifnot(all(u > 0 & u < 1))
    stopifnot(inherits(obj, "kdecopula"))
    stopifnot(is.logical(stable))
    
    ## define appropriately shaped udata frame for vapply
    u <- as.matrix(u)
    if (ncol(u) == 1)
        u <- matrix(u, 1L, nrow(u))
    
    ## adjust for flipping option of kdevine package
    d <- ncol(u)
    if (!is.null(obj$flip))
        u <- matrix(u[, 2:1], nrow(u), d)
    
    
    ## if independence copula is specified return 1
    if ("indep.copula" %in% class(obj))
        return(rep(1, nrow(u)))
    
    ## evaluate density  (use faster algorithm for d = 2)
    u <- pmin(pmax(u, 1e-3), 1 - 1e-3)
    if (d == 2) {
        out <- interp_2d(u,
                         obj$estimate,
                         obj$grid)
    } else {
        # define help indicators
        tmplst <- split(rep(seq(-1, 2, 1), d), ceiling(seq.int(4*d)/4))
        helpind <- as.matrix(do.call(expand.grid, tmplst))
        
        out <- interp(u,
                      obj$estimate,
                      obj$grid,
                      helpind)
    }
    
    ## stabilize output
    if (stable)
        out <- pmin(out, 10^(1 + d/2))
    
    ## return results
    out
}
