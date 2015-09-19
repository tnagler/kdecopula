
## plot kdecopula object


#' Plotting \code{kdecopula} objects
#' 
#' Produces contour or perspective plots for a \code{kdecopula} object.
#' 
#' 
#' @param x \code{kdecopula} object.
#' @param type plot type; either \code{"contour"} or \code{"surface"} (partial
#' matching is activated) for a contour or perspective/surface plot
#' respectively.
#' @param margins either \code{"unif"} for the original copula density or
#' \code{"norm"} for the transformed density with standard normal margins
#' (partial matching is activated). Default is \code{"norm"} for \code{type =
#' "contour"}, and \code{"unif"} for \code{type = "surface"}.
#' @param size integer; the plot is based on values on a \eqn{size x size}
#' grid; default is 100 for \code{type = "contour"}, and 25 for \code{type =
#' "surface"}.
#' @param \dots optional arguments passed to \code{\link{contour}} or
#' \code{\link{wireframe}}.
#' @author Thomas Nagler
#' @seealso \code{\link[kdecopula:kdecop]{kdecop}},
#' \code{\link[graphics:contour]{contour}},
#' \code{\link[lattice:wireframe]{wireframe}}
#' @keywords plot Copula density estimation
#' @examples
#' 
#' ## load data and transform with empirical cdf
#' data(wdbc)
#' udat <- apply(wdbc[, -1], 2, function(x) rank(x)/(length(x)+1))
#' 
#' ## estimation of copula density of variables 5 and 6
#' obj <- kdecop(udat[, 5:6])
#' 
#' ## plots
#' plot(obj)  # (marginal normal) contour plot
#' plot(obj, margins = "unif")  # contour plot of actual copula density
#' plot(obj, type = "surf")  # surface plot of actual copula density 
#' 
#' 
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
