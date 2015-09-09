kdecop <- function(udata, bw, mult = 1, method = "TLL2", knots = NA, renorm.iter = 3L, info = FALSE) {
    udata <- as.matrix(udata)
    n <- nrow(udata)
    d <- ncol(udata)
    
    ## sanity checks
    if (n < 2)
        stop("Number of observations has to be at least 2.")
    if (ncol(udata) < 2)
        stop("Dimension has to be at least 2.")
    if (any(udata > 1) | any(udata < 0))
        stop("'udata' have to be in the interval [0,1].")
    if (!(method %in% c("MR", "beta", "T", "TLL1", "TLL2")))
        stop("method not implemented")
    if (mult <= 0)
        stop("'mult' has to be a positive number.")
    if (is.na(knots))
        knots <- 100/d
    knots <- round(knots)
    stopifnot(is.numeric(renorm.iter))
    renorm.iter <- round(renorm.iter)
    stopifnot(renorm.iter >= 0)
    stopifnot(is.logical(info))
    
    ## bandwidth selection and adjustment (with bandwidth checks)
    if (missing(bw))
        bw <- bw_select(udata, method)
    if (method %in% c("TLL1", "TLL2")) {
        if (is.null(bw$B) | is.null(bw$alpha))
            stop("For methods 'TLL1/2', you have to provide a list 'bw = list(B = <your.B>,  alpha = <your.alpha>)'
where both parts of the bandwidth specification are provided via  'your.B', and 'your.alpha'.")
        if (any(c(diag(bw$B), bw$alpha) <= 0))
            stop("Bandwidths have to be positive.")
        bw$alpha <- mult * bw$alpha
    } else if (method %in% c("T")) {
        B <- as.matrix(bw)
        if (nrow(B) == 1)
            bw <- diag(as.numeric(bw), d)
        bw <- mult * bw
        bw <- bw * sqrt(min((n^(-1/6))^2 / det(bw), 1))
    } else if (method == "MR") {
        bw <- min(bw * mult, 1)
    } else if (method == "beta") {
        bw <- bw * mult
    }
    
    ## fit model for method TLL
    if (method %in% c("TLL1", "TLL2")) {
        zdata <-  qnorm(udata)
        lfit <- my_locfit(zdata,
                          bw$B,
                          bw$alpha,
                          deg = as.numeric(substr(method, 4, 4)))
    } else {
        lfit <- NA
    }
    
    ## construct grid with k knots in dimension d
    pnts <- pnorm(seq(-3.25, 3.25, l = knots))
    # pnts <- 1:knots/(knots + 1)
    grid <- as.matrix(do.call(expand.grid,
                              split(rep(pnts, d), ceiling(1:(knots*d)/knots))))
    
    ## evaluate estimator on grid
    evalf <- eval_func(method)  # get evaluation function
    object <- list(udata = udata,
                   bw = bw,
                   lfit = lfit,
                   method = method)
    vals <- array(evalf(grid, obj=object), dim = rep(knots, d))
    
    ## rescale copula density to have uniform margins
    if (renorm.iter > 0) {
        tmplst <- split(rep(seq.int(knots)-1, d-1),
                        ceiling(seq.int(knots*(d-1))/knots))
        helpind <- as.matrix(do.call(expand.grid, tmplst))
        vals <- renorm(vals, pnts, renorm.iter, helpind)
    }
    
    ## store results
    res <- list(udata    = udata,
                grid     = pnts,
                estimate = vals,
                bw       = bw,
                method   = method)
    class(res) <- "kdecopula"
    
    if (info) {
        # likelihood
        likvalues <- dkdecop(udata, res)
        loglik <- sum(log(likvalues))
        # effective number of parameters
        effp <- eff_num_par(udata, likvalues, bw, method, lfit)
        # information criteria
        AIC  <- - 2 * loglik + 2 * effp
        cAIC <- AIC + (2 * effp * (effp + 1)) / (n - effp - 1)
        BIC  <- - 2 * loglik + log(n) * effp
        
        ## store results
        res <- append(res, list(info = list(likvalues = likvalues,
                                            loglik    = loglik,
                                            effp      = effp ,
                                            AIC       = AIC,
                                            cAIC      = cAIC,
                                            BIC       = BIC)))
    }
    
    ## return results
    res
}

##### evaluate density of kdecopula object
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

##### evaluate h-function of kdecopula object
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

## function for preparing the list for constructing the helpgrid in hkdecop
prep_hfunc <- function(i, cond.var, grid, d) {
    if (i %in% cond.var) {
        return(0)
    } else {
        return(grid)
    }
}

##### evaluate cdf of kdecopula object
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

##### simulate from kdecopula object
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


##### inverse h-function
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



## plot kdecopula object
plot.kdecopula <- function(x, type = "contour", margins, size, ...) {
    if (length(dim(x$estimate)) > 2)
        stop("Plotting is only available for 2-dimensional densities.")
    
    ## partial matching and sanity check for type
    stopifnot(class(type) == "character")
    tpnms <- c("contour", "surface")
    type <- tpnms[pmatch(type, tpnms)]
    if (is.na(type))
        stop("type not implemented")
    
    ## choose margins if missing, else partial matching and sanity check
    if (missing(margins)) {
        margins <- switch(type,
                          "contour" = "norm",
                          "surface" = "unif")
    } else {
        stopifnot(class(margins) == "character")
        mgnms <- c("norm", "unif")
        margins <- mgnms[pmatch(margins, mgnms)]
    }
    
    ## choose size if missing and sanity check
    if (missing(size))
        size <- switch(type,
                       "contour" = 100L,
                       "surface" = 25L)
    stopifnot(is.numeric(size))
    
    
    ## construct grid for evaluation of the copula density
    size <- round(size)
    if (size < 3) {
        warning("size too small, set to 5")
        size <- 5
    }
    if (!(margins %in% c("unif", "norm")))
        stop("'margins' has to be one of 'unif' or 'norm'")
    if (is.null(list(...)$xlim) & is.null(list(...)$ylim)) {
        xylim <- switch(margins,
                        "unif"  = c(0, 1),
                        "norm"  = c(-3, 3))
    } else {
        xylim <- range(c(list(...)$xlim, list(...)$ylim))
    }
    sq <- seq(xylim[1L], xylim[2L], len = size)
    points <- switch(margins,
                     "unif"  = 1:size/(size + 1),
                     "norm"  = pnorm(sq))
    g <- as.matrix(expand.grid(points, points))
    
    ## evaluate on grid
    vals <- dkdecop(g, x)
    cop <- matrix(vals, size, size)
    
    ## prepare for plotting with selected margins
    if (margins == "unif") {
        points <- g[1L:size, 1L]
        adj <- 1
        gu <- g[, 1L]
        gv <- g[, 2L]
        levels <- c(0.1, 0.5, 1, 3, 5, 10, 20)
        xlim <- ylim <- c(0, 1)
        at <- c(seq(0, 6, by = 0.05), seq(7, 100, by = 1))
    } else if (margins == "norm") {
        points <- qnorm(g[1L:size, 1L])
        adj <- tcrossprod(dnorm(points))
        gu <- qnorm(g[, 1L])
        gv <- qnorm(g[, 2L])
        levels <- c(0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
        xlim <- ylim <- c(-3, 3)
        at <- seq(0, 1, l = 100)
    }
    
    if (type == "contour") {
        # set default parameters
        pars <- list(x = points,
                     y = points,
                     z = cop * adj,
                     levels = levels,
                     xlim = xlim,
                     ylim = ylim,
                     xlab = switch(margins,
                                   "unif" = expression(u[1]),
                                   "norm" = expression(z[1])),
                     ylab = switch(margins,
                                   "unif" = expression(u[2]),
                                   "norm" = expression(z[2])))
        
        # call contour with final parameters
        do.call(contour, modifyList(pars, list(...)))
        
    } else if (type == "heat") {
        stop("Not implemented yet")
    } else if (type == "surface") {
        # list with coordinates
        lst <- list(u = gu, v = gv, c = as.vector(cop) * as.vector(adj))
        
        # define colors
        TUMblue   <- rgb(0, 103/255, 198/255)
        TUMgreen  <- rgb(162/255, 173/255, 0)
        TUMorange <- rgb(227/255, 114/255, 37/255)
        
        # set default parameters
        pars <- list(x = c ~ u * v,
                     data = lst,
                     scales = list(arrows = FALSE),
                     drape = TRUE, colorkey = FALSE,
                     screen = list(z = 25, x = -55),
                     shade = FALSE,
                     aspect = c(1, 1),
                     light.source = c(10,0,10),
                     zoom = 0.85,
                     par.settings = list(axis.line = list(col = "transparent")),
                     at = at,
                     col.regions=
                         c(colorRampPalette(
                             c(TUMblue, TUMgreen, TUMorange))(121),
                           rep(TUMorange, 300)),
                     xlab = switch(margins,
                                   "unif" = expression(u[1]),
                                   "norm" = expression(z[1])),
                     ylab = switch(margins,
                                   "unif" = expression(u[2]),
                                   "norm" = expression(z[2])),
                     zlab = "density")
        
        # call wireframe with final parameters
        do.call(wireframe, modifyList(pars, list(...)))
    }
}
