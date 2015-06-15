### bandwidth selection
bw_select <- function(udata, method) {
    n <- nrow(udata)
    d <- ncol(udata)
    switch(method,
           "MR"   = bw_mr(udata),
           "beta" = bw_beta(udata),
           "T"    = n^(-1/(d + 4)) * t(chol(cov(qnorm(udata)))),
           "TLL1" = bw_tll(qnorm(udata), deg = 1),
           "TLL2" = bw_tll(qnorm(udata), deg = 2))
}

bw_mr <- function(udata) {
    n <- nrow(udata)
    
    ## constants for kernel
    sigma_K <- sqrt(1/5)
    d_K     <- 3/5
    
    ## parameter for frank copula by inversion of Kendall's tau
    family <- 5
    tau <- cor(udata, method="kendall")[1L, 2L]
    if (abs(tau) < 1e-16) {
        family <- 0
        par <- 0
    } else {
        par <- BiCopTau2Par(family, tau=tau)
    }
    
    ## short handles for copula density and derivatives
    cd   <- function(u,v) 
        BiCopPDF(u, v, family, par)
    c_u  <- function(u,v)
        BiCopDeriv(u, v, family, par, deriv="u1")
    c_v  <- function(u,v) 
        BiCopDeriv(u, v, family, par, deriv="u2")
    c_uu <- function(u,v) 
        BiCopDeriv2(u, v, family, par, deriv="u1")
    c_vv <- function(u,v) 
        BiCopDeriv2(u, v, family, par, deriv="u2")
    
    ## short handles for integrands
    bet <- function(w) (c_uu(w[1L], w[2L]) + c_vv(w[1L], w[2L]))^2
    gam <- function(w) cd(w[1L], w[2L])
    
    ## integrations
    beta  <- adaptIntegrate(bet,
                            lowerLimit = c(0,0),
                            upperLimit = c(1,1),
                            tol = 5e-3)$integral
    gamma <- adaptIntegrate(gam,
                            lowerLimit = c(0,0),
                            upperLimit = c(1,1),
                            tol = 5e-3)$integral
    
    ## result
    res <- (2*d_K^2/sigma_K^4*gamma/beta)^(1/6) * n^(-1/6)
    if (res > 1) 1 else res
}

bw_beta <- function(udata) {
    n  <- nrow(udata)
    
    ## parameter for frank copula by inversion of Kendall's tau
    family <- 5
    R <- cor(udata, method="kendall")
    tau <- mean(R[lower.tri(R)])
    if (abs(tau) < 1e-16) {
        family <- 0
        par <- 0
    } else {
        par <- BiCopTau2Par(family, tau=tau)
    }
    
    ## short handles for copula density and derivatives
    cd   <- function(u,v)
        BiCopPDF(u, v, family, par)
    c_u  <- function(u,v)
        BiCopDeriv(u, v, family, par, deriv = "u1")
    c_v  <- function(u,v)
        BiCopDeriv(u, v, family, par, deriv = "u2")
    c_uu <- function(u,v)
        BiCopDeriv2(u, v, family, par, deriv = "u1")
    c_vv <- function(u,v)
        BiCopDeriv2(u, v, family, par, deriv = "u2")
    
    ## short handles for integrands
    x <- function(w) {
        u <- w[1L]
        v <- w[2L]
        ((1-2*u)*c_u(u,v) + (1-2*v)*c_v(u,v) + 
             1/2 * (u*(1-u)*c_uu(u,v) + v*(1-v)*c_vv(u,v)))^2
    }
    zet <- function(w) {
        u <- w[1L]
        v <- w[2L]
        cd(u,v) / sqrt(u*(1-u)*v*(1-v))
    }
    
    ## integrations
    xi   <- adaptIntegrate(x,
                           lowerLimit = c(0,0),
                           upperLimit = c(1,1),
                           tol = 5e-3,
                           maxEval = 10^3)$integral
    zeta <- adaptIntegrate(zet,
                           lowerLimit = c(0,0), 
                           upperLimit = c(1,1),
                           tol = 5e-3,
                           maxEval = 10^3)$integral
    
    ## result
    (zeta/(8*pi*xi))^(1/3) * n^(-1/3)
}

bw_tll <- function(zdata, deg) {
    n <- nrow(zdata)
    d <- ncol(zdata)
    # transform to uncorrelated data
    B <- eigen(cov(zdata))$vectors
    B <- B * sign(diag(B))
    qrs <- zdata %*% B
    
    ## find optimal alpha in for each principal component
    # function to optimize over
    fn <- function(alpha, x) {
        e <- suppressWarnings(try(lscv(lp(x, nn = alpha, deg = deg),
                                       kern = "gauss"),
                                  silent = TRUE))
        if ("try-error" %in% class(e)) Inf else e[1]
    }
    # optimization function
    opt <- function(i) {
        optim(0.2, 
              fn, 
              x = qrs[, i],
              lower = 1e-6,
              upper = 1,
              method = "Brent")$par
    }
    alpha.vec <- sapply(1:d, opt)
    
    ## adjustments for multivariate estimation and transformation
    kappa <- alpha.vec[1]/alpha.vec
    B <- B %*% diag(1/kappa)
    dimnames(B) <- NULL
    if (deg == 1) {
        alpha <- n^(1/5 - d/(4 + d)) * alpha.vec[1]
    } else {
        alpha <- n^(1/9 - d/(8 + d)) * alpha.vec[1]
    }
    
    ## return results
    list(B = B, alpha = alpha)
}


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
        evalpoints <- cbind(rep(U/b, 9) - augdata[, 1]/b,
                            rep(V/b, 9) - augdata[, 2]/b)
        K <- kern_gauss_2d(evalpoints[, 1]/b, evalpoints[, 2]/b)
        S <- rowSums(matrix(K, n, 9)) / b^2
        effp <- sum(S/likvalues)/n
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
        effp  <- mean(kern_gauss_2d(0, 0)/scale/det(b)/likvalues)
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
           "TB"   = function(uev, obj) 
               eval_trafo(uev, obj$udata, obj$bw),
           "TLL1" = function(uev, obj) 
               eval_tll(uev, obj$lfit, obj$bw$B),
           "TLL2" = function(uev, obj) 
               eval_tll(uev, obj$lfit, obj$bw$B))
}

##### local likelihood fitting
my_locfit <- function(zdata, B, alpha, deg) {
    ## construct grid (first u level, then transform to z level with B)
    d <- ncol(zdata)
    m <- round(100/d)
    tmplst <- split(rep(seq.int(m)/(m+1), d), ceiling(seq.int(m*d)/m))
    gr   <- as.matrix(do.call(expand.grid, tmplst))
    grQR <- qnorm(gr) %*% B
    
    ## transform data
    qrs  <- zdata %*% B
    
    ## fit model
    lims <- apply(grQR, 2L, range) * 1.8
    cl.lst <- split(as.vector(qrs), ceiling(seq.int(nrow(qrs)*d)/nrow(qrs)))
    cl.lst$nn <- alpha
    cl.lst$deg <- deg
    lf.lst <- list(~do.call(lp, cl.lst),
                   maxk = 1000,
                   kern = "gauss")
    if (d == 2) {
        lf.lst$ev <- lfgrid(mg = m, 
                            ll = lims[1L, ], 
                            ur = lims[2L, ]) 
    }
    suppressWarnings(do.call(locfit, lf.lst))
}

