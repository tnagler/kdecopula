### effective number of parameters
eff_num_par <- function(udata, likvalues, b, method, lfit) {
    if (method %in% c("TLL1", "TLL2", "TLL1nn", "TLL2nn")) {
        effp <- lfit$dp[["df2"]]
    } else if (method %in% c("MR")) {
        U <- udata[, 1]
        V <- udata[, 2]
        n <- length(U)
        augdata <- rbind(cbind(U, V),
                         cbind(-U, V),
                         cbind(U, -V),
                         cbind(-U, -V),
                         cbind(U, 2 - V),
                         cbind(-U, 2 - V),
                         cbind(2 - U, V),
                         cbind(2 - U, -V),
                         cbind(2 - U, 2 - V))
        evalpoints <- cbind(rep(U, 9) - augdata[, 1],
                            rep(V, 9) - augdata[, 2])
        K <- kern_gauss_2d(evalpoints[, 1], evalpoints[, 2], b)
        S <- rowSums(matrix(K, n, 9))
        effp <- sum(S / likvalues) / n
    } else if (method == "beta") {
        bpkern <- function(x) {
            dbeta(x[1L], x[1L]/b + 1, (1-x[1L])/b + 1) *
                dbeta(x[2L], x[2L]/b + 1, (1-x[2L])/b + 1)
        }
        p <- apply(udata, 1, bpkern)
        effp <- mean(p/likvalues)
    } else if (method == "T") {
        scale <- dnorm(qnorm(udata)[, 1]) * dnorm(qnorm(udata)[, 2])
        effp  <- mean((kern_gauss_2d(0, 0, 1) / (scale * det(b))) / likvalues)
    }
    
    ## return result
    effp
}

##### return functions for evaluators
eval_func <- function(method) {
    switch(method,
           "MR"   = function(uev, obj)
               eval_mr(uev, obj$udata, obj$bw),
           "beta" = function(uev, obj)
               eval_beta(uev, obj$udata, obj$bw),
           "T"    = function(uev, obj)
               eval_trafo(uev, obj$udata, obj$bw),
           "TLL1" = function(uev, obj)
               eval_tll(uev, obj$lfit, obj$bw),
           "TLL2" = function(uev, obj)
               eval_tll(uev, obj$lfit, obj$bw),
           "TLL1nn" = function(uev, obj)
               eval_tll(uev, obj$lfit, obj$bw$B),
           "TLL2nn" = function(uev, obj)
               eval_tll(uev, obj$lfit, obj$bw$B),
           "TTPI" = function(uev, obj)
               eval_tt(uev, obj$udata, obj$bw),
           "TTCV" = function(uev, obj)
               eval_tt(uev, obj$udata, obj$bw),
           "bern" = function(uev, obj)
               dberncop(uev, berncop(obj$udata, obj$bw)))
}



#' Title
#'
#' @param udata 
#' @param B 
#' @param deg 
#'
#' @return
#' 
#' @import locfit
my_locfit <- function(udata, B, deg) {
    if (is.list(B))
        return(my_locfitnn(udata, B$B, B$alpha, B$kappa, deg))
    # transform data
    zdata <- qnorm(udata)
    qrs  <- zdata %*% solve(B)
    gr <- do.call(expand.grid, lapply(1:ncol(qrs), function(i) c(-4, 4)))
    qgr <- as.matrix(gr) %*% solve(B)
    lims <- apply(qgr, 2L, range)
    
    ## fit model
    cl.lst <- split(as.vector(qrs), rep(1:ncol(qrs), each = nrow(qrs)))
    cl.lst$h <- 1
    cl.lst$deg <- deg
    lf.lst <- list(~do.call(lp, cl.lst),
                   maxk = 1000,
                   kern = "gauss",
                   scale = FALSE,
                   ev = lfgrid(mg = 50, 
                               ll = lims[1L, ], 
                               ur = lims[2L, ]))
    suppressWarnings(do.call(locfit, lf.lst))
}

my_locfitnn <- function(udata, B, alpha, kappa, deg) {
    # transform data
    zdata <- qnorm(udata)
    qrs  <- zdata %*% solve(B)
    gr <- do.call(expand.grid, lapply(1:ncol(qrs), function(i) c(-4, 4)))
    qgr <- as.matrix(gr) %*% solve(B)
    lims <- apply(qgr, 2L, range)
    
    ## fit model
    cl.lst <- split(as.vector(qrs), rep(1:ncol(qrs), each = nrow(qrs)))
    cl.lst$nn <- alpha
    cl.lst$deg <- deg
    cl.lst$scale <- kappa
    lf.lst <- list(~do.call(lp, cl.lst),
                   maxk = 1024,
                   kern = "gauss",
                   ev = lfgrid(mg = 50, 
                               ll = lims[1L, ], 
                               ur = lims[2L, ]))
    suppressWarnings(do.call(locfit, lf.lst))
}


