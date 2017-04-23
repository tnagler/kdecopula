#' @export
print.kdecopula <- function(x, ...) {
    cat("Kernel copula density estimate '", x$method, "'", sep = "")
    ## add variable names if available
    nms <- colnames(x$udata)
    if (length(nms) == 2) {
        cat(":", nms[1], "--", nms[2])
    }
    invisible(x)
}



#' @importFrom stats logLik
#' @export
summary.kdecopula <- function(object, ...) {
    cat("Kernel copula density estimate\n")
    cat("------------------------------\n")

    ## add variable names if available
    nms <- colnames(object$udata)
    if (length(nms) == 2) {
        cat("Variables:   ", nms[1], "--", nms[2], "\n")
    }

    ## more details about the estimate
    cat("Observations:",
        nrow(object$udata), "\n")
    cat("Method:      ",
        expand_method(object$method), "\n")
    if (object$method %in% c("MR", "beta")) {
        cat("Bandwidth:   ",
            round(object$bw,3), "\n")
    } else if (object$method %in% c("TLL1nn", "TLL2nn")) {
        cat("Bandwidth:    ",
            "alpha = ",
            object$bw$alpha, "\n              ",
            "B = matrix(c(",
            paste(round(object$bw$B, 2), collapse = ", "),
            "), 2, 2)\n",
            sep = "")
    } else if (object$method %in% c("TLL1", "TLL2")) {
        cat("Bandwidth:    ",
            "matrix(c(",
            paste(round(object$bw, 2), collapse = ", "),
            "), 2, 2)\n",
            sep = "")
    } else if (object$method %in% c("T")) {
        cat("Bandwidth:    ",
            "matrix(c(", paste(round(object$bw, 2), collapse = ", "), "), 2, 2)\n",
            sep = "")
    } else if (object$method %in% c("TTPI", "TTCV")) {
        cat("Bandwidth:    ",
            "(",
            paste(round(object$bw, 2), collapse = ","),
            ")\n",
            sep = "")
    }
    cat("---\n")

    ## fit statistics
    if (object$method %in% c("TTPI", "TTCV")) {
        cat("logLik:", round(logLik(object), 2), "\n")
        cat("No further information available for this method")
    } else {
        # calculate
        effp <- attr(logLik(object), "df")
        loglik <- logLik(object)
        AIC  <- AIC(object)
        cAIC <- AIC + (2 * effp * (effp + 1)) / (nrow(object$udata) - effp - 1)
        BIC  <- BIC(object)
        # store in object if not available
        if (is.null(object$info))  object$info <- list(loglik = as.numeric(loglik),
                                                       effp   = effp ,
                                                       AIC    = AIC,
                                                       cAIC   = cAIC,
                                                       BIC    = BIC)
        # print
        cat("logLik:",
            round(loglik, 2), "   ")
        cat("AIC:",
            round(AIC, 2), "   ")
        cat("cAIC:",
            round(cAIC, 2), "   ")
        cat("BIC:",
            round(BIC, 2), "\n")
        cat("Effective number of parameters:",
            round(effp, 2))
    }
    ## return with fit statistics
    invisible(object)
}

expand_method <- function(method) {
    switch(method,
           "TLL1"   = "Transformation local likelihood, log-linear ('TLL1')",
           "TLL2"   = "Transformation local likelihood, log-quadratic ('TLL2')",
           "TLL1nn" = "Transformation local likelihood, log-linear (nearest-neighbor, 'TLL1nn')",
           "TLL2nn" = "Transformation local likelihood, log-quadratic (nearest-neighbor, 'TLL2nn')",
           "T"      = "Transformation estimator ('T')",
           "MR"     = "Mirror-reflection ('MR')",
           "beta"   = "Beta kernels ('beta')",
           "TTPI"   = "Tapered transformation estimator (plug-in, 'TTPI')",
           "TTCV"   = "Tapered transformation estimator (cross-validated, 'TTCV')")
}


#' Log-Likelihood of a \code{kdecopula} object
#'
#' @method logLik kdecopula
#'
#' @param object an object of class \code{kdecopula}.
#' @param ... not used.
#'
#' @return Returns an object of class \code{\link[stats:logLik]{logLik}} containing the log-
#' likelihood, number of observations and effective number of parameters ("df").
#'
#' @author Thomas Nagler
#'
#' @seealso
#' \code{\link[stats:logLik]{logLik}},
#' \code{\link[stats:AIC]{AIC}},
#' \code{\link[stats:BIC]{BIC}}
#'
#' @examples
#' ## load data and transform with empirical cdf
#' data(wdbc)
#' udat <- apply(wdbc[, -1], 2, function(x) rank(x)/(length(x)+1))
#'
#' ## estimation of copula density of variables 5 and 6
#' dens.est <- kdecop(udat[, 5:6])
#'
#' ## compute fit statistics
#' logLik(dens.est)
#' AIC(dens.est)
#' BIC(dens.est)
#'
#' @export
logLik.kdecopula <- function(object, ...) {
    if (!is.null(object$info)) {
        ## access info slot if available
        out <- object$info$loglik
        effp <- object$info$effp
    } else {
        ## calculate log likelihood and effp from data
        likvalues <- dkdecop(object$udata, object)
        out <- sum(log(likvalues))
        effp <- eff_num_par(object$udata,
                            likvalues,
                            object$bw,
                            object$method,
                            object$lfit)
    }

    ## add attributes
    attr(out, "nobs") <- nrow(object$udata)
    attr(out, "df") <- effp

    ## return object with class "logLik"
    class(out) <- "logLik"
    out
}

#' Prediction method for `kdecop()` fits
#' 
#' Predicts the pdf, cdf, or (inverse) h-functions by calling `dkdecop()`, 
#' `pkdecop()`, or `hkdecop()`.
#'
#' @param object an object of class `kdecopula`.
#' @param newdata evaluation point(s), a length two vector or a matrix with
#'   two columns.
#' @param what what to predict, one of `c("pdf", "cdf", "hfunc1", "hfunc2",
#' "hinv1", "hinv2")`. 
#' @param stable only used for `what = "pdf"`, see `dkdecop()`.
#' @param ... unused.
#'
#' @return A numeric vector of predictions.
#' 
#' @examples
#' data(wdbc)
#' udat <- apply(wdbc[, -1], 2, function(x) rank(x)/(length(x)+1))
#' est <- kdecop(udat[, 5:6])
#' 
#' all.equal(predict(est, c(0.1, 0.2)), dkdecop(c(0.1, 0.2), est))
#' all.equal(predict(est, udat, "hfunc1"), hkdecop(udat, est, cond.var = 1))
#' @export
#' @importFrom stats predict
predict.kdecopula <- function(object, newdata, what = "pdf", stable = FALSE, ...) {
    switch(
        what,
        "pdf"    = dkdecop(newdata, object, stable),
        "cdf"    = pkdecop(newdata, object),
        "hfunc1" = hkdecop(newdata, object, cond.var = 1),
        "hfunc2" = hkdecop(newdata, object, cond.var = 2),
        "hinv1"  = hkdecop(newdata, object, cond.var = 1, inverse = TRUE),
        "hinv2"  = hkdecop(newdata, object, cond.var = 2, inverse = TRUE)
    )
}

#' Extract fitted values from a `kdecop()` fits.
#' 
#' Simply calls `predict(object, object$udata, what)`.
#'
#' @param object an object of class `kdecopula`.
#' @param what what to predict, one of `c("pdf", "cdf", "hfunc1", "hfunc2",
#' "hinv1", "hinv2")`. 
#' @param ... unused.
#'
#' @seealso `predict.kdecopula()`
#' @examples
#' data(wdbc)
#' udat <- apply(wdbc[, -1], 2, function(x) rank(x)/(length(x)+1))
#' est <- kdecop(udat[, 5:6])
#' 
#' all.equal(fitted(est), predict(est, est$udata))
#' @export
#' @importFrom stats fitted
fitted.kdecopula <- function(object, what = "pdf", ...) {
    predict(object, object$udata, what)
}

#' Simulate synthetic data from a `kdecop()` fits.
#' 
#' See `rkdecop()`.
#' 
#'
#' @param object an object of class `kdecopula`.
#' @param nsim integer; number of observations.
#' @param seed integer; `set.seed(seed)` will be called prior to `rkdecop()`.
#' @param quasi logical; the default (\code{FALSE}) returns pseudo-random
#' numbers, use \code{TRUE} for quasi-random numbers (generalized Halton, see
#' @param ... unused.
#'
#' @return Simulated data from the fitted `kdecopula` model.
#'
#' @examples
#' 
#' @export 
simulate.kdecopula <- function(object, nsim = 1, seed = NULL, quasi = FALSE, ...) {
    if (!is.null(seed))
        set.seed(seed)
    rkdecop(nsim, object, quasi)
}



