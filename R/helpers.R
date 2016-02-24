### effective number of parameters
eff_num_par <- function(udata, likvalues, b, method, lfit) {
    if (method %in% c("MR")) {
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
    }
    if (method == "beta") {
        bpkern <- function(x) {
            dbeta(x[1L], x[1L]/b + 1, (1-x[1L])/b + 1) *
                dbeta(x[2L], x[2L]/b + 1, (1-x[2L])/b + 1)
        }
        p <- apply(udata, 1, bpkern)
        effp <- mean(p/likvalues)
    }
    if (method == "T") {
        scale <- dnorm(qnorm(udata)[, 1]) * dnorm(qnorm(udata)[, 2])
        effp  <- mean((kern_gauss_2d(0, 0, 1) / (scale * det(b))) / likvalues)
    }
    if(method %in% c("TLL1", "TLL2"))
        effp <- lfit$dp[["df2"]]

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
               eval_tll(uev, obj$lfit, obj$bw$B),
           "TLL2" = function(uev, obj)
               eval_tll(uev, obj$lfit, obj$bw$B),
           "TTPI" = function(uev, obj)
               eval_tt(uev, obj$udata, obj$bw),
           "TTCV" = function(uev, obj)
               eval_tt(uev, obj$udata, obj$bw))
}

##### local likelihood fitting
my_locfit <- function(zdata, B, alpha, deg) {
    ## construct grid (first u level, then transform to z level with B)
    d <- ncol(zdata)
    #     m <- round(100/d)
    #     tmplst <- split(rep(seq.int(m)/(m+1), d), ceiling(seq.int(m*d)/m))
    #     gr   <- as.matrix(do.call(expand.grid, tmplst))
    #     grQR <- qnorm(gr) %*% B

    ## transform data
    qrs  <- zdata %*% B

    ## fit model
    #     lims <- apply(grQR, 2L, range) * 1.8
    cl.lst <- split(as.vector(qrs), ceiling(seq.int(nrow(qrs)*d)/nrow(qrs)))
    cl.lst$nn <- alpha
    cl.lst$deg <- deg
    lf.lst <- list(~do.call(lp, cl.lst),
                   maxk = 1000,
                   kern = "gauss")
    suppressWarnings(do.call(locfit, lf.lst))
}

