#' Bivariate kernel copula density estimation
#' 
#' Based on samples from a bivariate copula, the copula density is estimated.
#' The user can choose between different methods. If no bandwidth is provided
#' by the user, it will be set by a method-specific automatic selection
#' procedure. The related (d/p/r)kdecop functions evaluate the density and cdf 
#' or simulate synthetic data, respectively.
#'   
#' @param udata \code{nx2} matrix of copula data.
#' @param bw bandwidth specification; if \code{NA}, \code{bw} is selected
#' automatically (see Details); Otherwise, please provide \cr
#' \code{"T", "TLL1", "TLL2"}: a \eqn{2x2} bandwidth matrix, \cr
#' \code{"TLL1nn", "TLL2nn"}: a list with (named) entries \code{B}, \code{alpha},
#' and \code{kappa}, \cr
#' \code{"TTCV", "TTPI"}: a numeric vector of length four containing \eqn{(h, 
#' \rho, \theta_1, \theta_2)}, c.f. Wen and Wu (2015), \cr
#' \code{"MR", "beta"}: a positive real number.
#' @param mult bandwidth multiplier, has to be positive; useful for making 
#' estimates more/less smooth manually.
#' @param method 
#' \code{"T"}: transformation estimator based on classical bivariate kernel 
#' estimation (e.g., Geenenens et al., 2014), \cr 
#' \code{"TLL1"}: transformation estimator with log-linear local likelihood
#' estimation (Geenenens et al., 2014), \cr 
#' \code{"TLL2"}: transformation estimator with log-quadradtic local likelihood
#' estimation (Geenenens et al., 2014), \cr 
#' \code{"TLL1nn"}: transformation estimator with log-linear local likelihood
#' estimation and nearest-neighbor bandwidths (Geenenens et al., 2014), \cr 
#' \code{"TLL2nn"}: transformation estimator with log-quadradtic local likelihood 
#' estimation and nearest-neighbor bandwidths (Geenenens et al., 2014), \cr
#' \code{"TTPI"}: tapered transformation estimator with plug-in bandwidths
#' (Wu and Wen, 2015), \cr
#' \code{"TTCV"}: tapered transformation estimator with profile cross-validation
#' bandwidths (Wu and Wen, 2015), \cr
#' \code{"MR"}: mirror-reflection estimator (Gijbels and Mielniczuk, 1990), \cr 
#' \code{"beta"}: beta kernel estimator (Charpentier et al., 2006), \cr
#' \code{"bern"}: Bernstein copula estimator (Sanchetta and Satchell, 2004); the
#' coefficients are adjusted by the method of Weiss and Scheffer (2012).
#' @param knots integer; number of knots in each dimension for the spline
#' approximation.
#' @param renorm.iter integer; number of iterations for the renormalization
#' procedure (see \emph{Details}).
#' @param info logical; if \code{TRUE}, additional information about the
#' estimate will be gathered (see \emph{Value}).
#' 
#' 
#' @return The function \code{\link[kdecopula:kdecop]{kdecop}} returns an
#' object of class \code{kdecopula} that contains all information necessary for
#' evaluation of the estimator. If no bandwidth was provided in the function
#' call, the automatically selected value can be found in the variable
#' \code{object$bw}. If \code{info=TRUE}, also the following will be available
#' under \code{object$info}: 
#' \item{likvalues}{Estimator evaluated in sample points} 
#' \item{loglik}{Log likelihood} 
#' \item{effp}{Effective number of parameters} 
#' \item{AIC}{Akaike information criterion}
#' \item{cAIC}{Bias-corrected version of Akaike information criterion}
#' \item{BIC}{Bayesian information criterion.} \cr 
#' The density estimate can be evaluated on arbitrary points with 
#' \code{\link[kdecopula:dkdecop]{dkdecop}}; the cdf with 
#' \code{\link[kdecopula:pkdecop]{pkdecop}}. Furthermore, synthetic data can be
#' simulated with \code{\link[kdecopula:rkdecop]{rkdecop}}, and several plotting
#' options are available with \code{\link[kdecopula:plot.kdecopula]{plot}}
#' and \code{\link[kdecopula:contour.kdecopula]{contour}}.
#' 
#' @details We use a Gaussian product kernel function for all methods 
#' except the beta kernel and Bernstein estimators. For details on bandwidth 
#' selection for a specific method, see: \code{\link{bw_t}}, 
#' \code{\link{bw_tll}}, \code{\link{bw_tll_nn}}, \code{\link{bw_tt_pi}}, 
#' \code{\link{bw_tt_cv}}, \code{\link{bw_mr}}, \code{\link{bw_beta}}, 
#' \code{\link{bw_bern}}.
#' \cr 
#' 
#' Kernel estimates are usually no proper copula densities. In particular, the
#' estimated marginal densities are not uniform. We mitigate this issue by
#' a renormalization procedure. The number of iterations of the
#' renormalization algorithm can be specified with the \code{renorm.iter}
#' argument. Typically, a very small number of iterations is sufficient. \cr
#' 
#' @note 
#' The implementation of the tapered transformation estimator ("TTPI"/"TTCV") 
#' was kindly provided by Kuangyu Wen. 
#' 
#' @author Thomas Nagler
#' 
#' @seealso 
#' \code{\link[kdecopula:kdecopula]{kdecopula}},
#' \code{\link[kdecopula:plot.kdecopula]{plot.kdecopula}},
#' \code{\link[kdecopula:dkdecop]{dkdecop}},
#' \code{\link[kdecopula:pkdecop]{pkdecop}},
#' \code{\link[kdecopula:rkdecop]{rkdecop}}
#' 
#' @references 
#' Geenens, G., Charpentier, A., and Paindaveine, D. (2017).
#' Probit transformation for nonparametric kernel estimation of the copula
#' density.
#' Bernoulli, 23(3), 1848-1873. 
#' \cr \cr
#' Wen, K. and Wu, X. (2015).
#' Transformation-Kernel Estimation of the Copula Density,
#' Working paper,
#' \url{http://agecon2.tamu.edu/people/faculty/wu-ximing/agecon2/public/copula.pdf}
#' \cr \cr
#' Gijbels, I. and Mielniczuk, J. (1990).
#' Estimating the density of a copula function.
#' Communications in Statistics - Theory and Methods, 19(2):445-464. 
#' \cr \cr 
#' Charpentier, A., Fermanian, J.-D., and Scaillet, O. (2006).
#' The estimation of copulas: Theory and practice. 
#' In Rank, J., editor, Copulas: From theory to application in finance. Risk Books.
#' \cr \cr 
#' Weiss, G. and Scheffer, M. (2012).
#' Smooth Nonparametric Bernstein Vine Copulas.
#' arXiv:1210.2043 [q-fin.RM]
#' \cr \cr
#' Nagler, T. (2014). 
#' Kernel Methods for Vine Copula Estimation.
#' Master's Thesis, Technische Universitaet Muenchen,
#' \url{https://mediatum.ub.tum.de/node?id=1231221} 
#' 
#' @examples
#' 
#' ## load data and transform with empirical cdf
#' data(wdbc)
#' udat <- apply(wdbc[, -1], 2, function(x) rank(x)/(length(x)+1))
#' 
#' ## estimation of copula density of variables 5 and 6
#' est <- kdecop(udat[, 5:6])
#' summary(est)
#' plot(est) 
#' contour(est)
#' 
#' ## evaluate density estimate at (u1,u2)=(0.123,0.321)
#' dkdecop(c(0.123, 0.321), est) 
#' 
#' ## evaluate cdf estimate at (u1,u2)=(0.123,0.321)
#' pkdecop(c(0.123, 0.321), est) 
#' 
#' ## simulate 500 samples from density estimate
#' plot(rkdecop(500, est))  # pseudo-random
#' plot(rkdecop(500, est), quasi = TRUE)  # quasi-random
#' 
#' 
#' @import lattice 
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' @importFrom graphics contour plot
#' 
#' @export
#' 
kdecop <- function(udata, bw = NA, mult = 1, method = "TLL2", knots = 30, 
                   renorm.iter = 3L, info = TRUE) {
    udata <- as.matrix(udata)
    n <- NROW(udata)
    d <- NCOL(udata)
    check_kdecop_input(as.list(environment()))
    
    # bandwidth selection
    if (any(is.na(bw)))
        bw <- bw_select(udata, method)
    bw <- multiply_bw(bw, mult, method, d)
    check_bw(bw, method)
    
    # fit model for method TLL
    if (method %in% c("TLL1", "TLL2", "TLL1nn", "TLL2nn")) {
        lfit <- my_locfit(udata, bw, deg = as.numeric(substr(method, 4, 4)))
    } else {
        lfit <- NULL
    }
    
    # evaluate estimator on grid
    grid <- make_grid(knots, d)
    evalf <- eval_func(method)
    object <- list(udata = udata,
                   bw = bw,
                   lfit = lfit,
                   method = method)
    vals <- array(evalf(grid$expanded, obj = object), dim = rep(knots, d))
    rm("object")
    # rescale copula density to have uniform margins
    vals <- renorm2unif(vals, grid, renorm.iter)
    
    # store results
    res <- list(udata    = udata,
                grid     = grid$pnts,
                estimate = vals,
                bw       = bw,
                mult     = mult,
                method   = method)
    class(res) <- "kdecopula"
    # add further (information if asked for)
    with_fit_info(res, info, lfit)
}

#' Sanity checks for kdecop function
#'
#' @param args list of all argument-value pairs
#'
#' @noRd
check_kdecop_input <- function(args) {
    if (args$n < 2)
        stop("Number of observations has to be at least 2.")
    if (args$d != 2)
        stop("Dimension has to be 2.")
    if (any(args$udata > 1) | any(args$udata < 0))
        stop("'udata' have to be in the interval [0,1].")
    
    if (!(args$method %in% c("MR", "beta", "T",
                             "TLL1", "TLL2", "TLL1nn", "TLL2nn",
                             "TTPI", "TTCV", "bern")))
        stop("method not implemented")
    
    if (args$mult <= 0)
        stop("'mult' has to be a positive number.")
    
    stopifnot(is.numeric(args$knots))
    if (!all.equal(args$knots, as.integer(args$knots)))
        stop("'knots' must be a positive integer.")
    if (args$knots < 1)
        stop("'knots' must be a positive integer.")
    
    stopifnot(is.numeric(args$renorm.iter))
    if (!all.equal(args$renorm.iter, as.integer(args$renorm.iter)))
        stop("'renorm.iter' must be a positive integer.")
    if (args$renorm.iter < 0)
        stop("'renorm.iter' must be a non-negative integer.")
    
    stopifnot(is.logical(args$info))
}

#' Grid for spline representation of the copula density estimate
#'
#' @param knots number of knots
#' @param d dimension
#' @noRd
make_grid <- function(knots, d) {
    pnts <- pnorm(seq(-3.25, 3.25, l = knots))
    pntlst <- split(rep(pnts, d), rep(c(1, 2), each = knots))
    expanded <- as.matrix(do.call(expand.grid, pntlst))
    list(pnts = pnts, expanded = expanded)
}

#' Normalizing the estimate to uniform margins
#'
#' @param vals initial values of the copula density evaluated on the grid
#' @param grid a list with expanded grid and univarite gridpoints
#' @param renorm.iter integer; number of iterations for the renormalization
#' procedure
#'
#' @return normalized values on the grid
#'
#' @noRd
renorm2unif <- function(vals, grid, renorm.iter) {
    if (renorm.iter > 0) {
        d <- NCOL(grid$expanded)
        knots <- length(grid$pnts)
        tmplst <- split(rep(seq.int(knots)-1, d-1),
                        ceiling(seq.int(knots*(d-1))/knots))
        helpind <- as.matrix(do.call(expand.grid, tmplst))
        vals <- renorm(vals, grid$pnts, renorm.iter, helpind)
    }
    vals
}

#' Finalize kdecopula object
#'
#' @param res kdecopula object
#' @param info logical; whether further information about the fit shall be 
#' computed
#' @param lfit locfit object (only for TLL methods, otherwise NULL)
#'
#' @return the kdecopula object with added info entry (if demanded).
#'
#' @noRd
with_fit_info <- function(res, info, lfit) {
    if (info) {
        # log-likelihood
        likvalues <- dkdecop(res$udata, res)
        loglik <- sum(log(likvalues))
        if (!(res$method %in% c("TTPI", "TTCV", "bern"))) {
            # effective number of parameters
            effp <- eff_num_par(res$udata,
                                likvalues, 
                                res$bw,
                                res$method, 
                                lfit)
            # information criteria
            AIC  <- - 2 * loglik + 2 * effp
            n <- nrow(res$udata)
            cAIC <- AIC + (2 * effp * (effp + 1)) / (n - effp - 1)
            BIC  <- - 2 * loglik + log(n) * effp
        } else {
            effp <- AIC <- cAIC <- BIC <- NA
        }
        
        ## store results
        res$info <- list(likvalues = likvalues,
                         loglik    = loglik,
                         effp      = effp ,
                         AIC       = AIC,
                         cAIC      = cAIC,
                         BIC       = BIC)
    }
    
    res
}
