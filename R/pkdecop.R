#' @rdname kdecop
#' 
pkdecop <- function(u, obj) {
    stopifnot(all(u > 0 & u < 1))
    stopifnot(inherits(obj, "kdecopula"))
    ## define appropriately shaped u matrix
    u <- as.matrix(u)
    if(ncol(u) == 1)
        u <- matrix(u, 1L, nrow(u))
    # adjust for flipping option of kdevine package
    if(!is.null(obj$flip)) {
        u <- matrix(u[, 2:1], nrow(u))
        cond.var <- ifelse(cond.var == 1, 2, 1)
    }
    d <- ncol(u)
    
    ## if independence copula is specified, return prod(u) directly
    if ("indep.copula" %in% class(obj))
        return(apply(u, 1, prod))
    
    ## define help objects
    tmplst <- split(rep(seq(-1, 2, 1), d), ceiling(seq.int(4*d)/4))
    helpind <- as.matrix(do.call(expand.grid, tmplst))
    m <- length(obj$grid)
    tmplst <- split(rep(obj$grid, d), ceiling(seq.int(m*d)/m))
    helpgrid <- as.matrix(do.call(expand.grid, tmplst))
    
    ## evaluate cdf
    eval_cdf(u,
             obj$estimate,
             obj$grid,
             helpgrid,
             helpind)
}
