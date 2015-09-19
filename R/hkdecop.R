# Evaluate h-function of kdecopula object (internal for kdevine-package)
# 
# @param cond.var integer; either 1 or 2; specifies on which argument to
#  condition 
# 
hkdecop <- function(u, obj, cond.var) {
    stopifnot(all(u > 0) & all(u < 1))
    stopifnot(inherits(obj, "kdecopula"))
    stopifnot(all(cond.var %in% c(1, 2)))
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
    
    ## if independence copula is specified, return the conditioned variable
    if ("indep.copula" %in% class(obj))
        return(u[, -cond.var])
    
    ## evaluate h-function otherwise
    if (d == 2) {
        # faster algorithm for d = 2
        return(eval_hfunc_2d(u,
                             as.integer(cond.var),
                             obj$estimate,
                             obj$grid))
    } else {
        # define help objects
        tmplst <- split(rep(seq(-1, 2, 1), d), ceiling(seq.int(4*d)/4))
        helpind <- as.matrix(do.call(expand.grid, tmplst))
        tmplst <- lapply(1:d,
                         prep_hfunc,
                         cond.var = cond.var,
                         grid = obj$grid,
                         d = d)
        helpgrid <- as.matrix(do.call(expand.grid, tmplst))
        uncond.var <- seq.int(d)[-cond.var]
        
        # call routine
        return(eval_hfunc(u,
                          as.integer(cond.var),
                          as.integer(uncond.var),
                          obj$estimate,
                          obj$grid,
                          helpgrid,
                          helpind))
    }
}

## Prepare the list for constructing the helpgrid in hkdecop
prep_hfunc <- function(i, cond.var, grid, d) {
    if (i %in% cond.var) {
        return(0)
    } else {
        return(grid)
    }
}

## Calculate inverse of h-function
hinvkdecop <- function(u, obj, cond.var) {
    stopifnot(all(u > 0 & u < 1))
    stopifnot(inherits(obj, "kdecopula"))
    stopifnot(all(cond.var %in% c(1, 2)))
    ## simulate independent uniform random variables
    W <- u
    
    # if independence copula is specified, return W
    if("indep.copula" %in% class(obj))
        return(W)
    
    ## invert h-function otherwise
    U2 <- inv_hfunc(W,
                    cond.var,
                    obj$estimate,
                    obj$grid)
    
    ## return results
    out <- cbind(W[, 1], U2)
    colnames(out) <- NULL
    out
}
