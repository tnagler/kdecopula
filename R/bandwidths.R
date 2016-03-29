### bandwidth selection
bw_select <- function(udata, method) {
    
    switch(method,
           "MR"   = bw_mr(udata),
           "beta" = bw_beta(udata),
           "T"    = nrow(udata)^(-1/(ncol(udata) + 4)) * t(chol(cov(qnorm(udata)))),
           "TLL1" = bw_tll(qnorm(udata), deg = 1),
           "TLL2" = bw_tll(qnorm(udata), deg = 2),
           "TTPI" = bw_tt_plugin(udata),
           "TTCV" = bw_tt_pcv(udata))
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
        par <- BiCopTau2Par(family, tau = tau)
    }
    
    ## short handles for copula density and derivatives
    c_uu <- function(u,v)
        BiCopDeriv2(u, v, family, par, deriv = "u1")
    c_vv <- function(u,v)
        BiCopDeriv2(u, v, family, par, deriv = "u2")
    
    ## integrals
    bet <- function(w) (c_uu(w[1L], w[2L]) + c_vv(w[1L], w[2L]))^2
    beta  <- adaptIntegrate(bet,
                            lowerLimit = c(0, 0),
                            upperLimit = c(1, 1),
                            tol = 5e-3)$integral
    gamma <- 1
    
    ## result
    res <- (2*d_K^2/sigma_K^4*gamma/beta)^(1/6) * n^(-1/6)
    if (res > 1) 1 else res
}

bw_beta <- function(udata) {
    n  <- nrow(udata)
    
    ## parameter for frank copula by inversion of Kendall's tau
    family <- 5
    R <- cor(udata, method = "kendall")
    tau <- mean(R[lower.tri(R)])
    if (abs(tau) < 1e-16) {
        family <- 0
        par <- 0
    } else {
        par <- BiCopTau2Par(family, tau = tau)
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
                           lowerLimit = c(0, 0),
                           upperLimit = c(1, 1),
                           tol = 5e-3,
                           maxEval = 10^3)$integral
    zeta <- adaptIntegrate(zet,
                           lowerLimit = c(0, 0),
                           upperLimit = c(1, 1),
                           tol = 5e-3,
                           maxEval = 10^3)$integral
    
    ## result
    (zeta/(8*pi*xi))^(1/3) * n^(-1/3)
}

bw_t <- function(udata) {
    ## normal rederence rule
    n <- nrow(udata)
    d <- ncol(udata)
    n^(- 1 / (d + 4)) * t(chol(cov(qnorm(udata)))) * (4 / (d + 2))^(1 / (d + 4))
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
        if (inherits(e, "try-error")) Inf else e[1]
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

## tapered transformation estimator
bw_tt_plugin <- function(obs, rho.add = T) {
    # This function uses the plug in method to select
    # the optimal smoothing parameters.  rho.add = T
    # indicates using the bandwidth matrix H = h^2 *
    # h^2 * (1, rho \\ rho, 1).  rho.add = F
    # indicates using the bandwidth matrix H= h^2 * (1,
    # 0 \\ 0, 1), namely the product kernel.
    n <- dim(obs)[1]
    Si <- qnorm(obs[, 1])
    Ti <- qnorm(obs[, 2])
    
    pre_rho <- mean(Si * Ti)
    b <- (32 * (1 - pre_rho^2)^3 * sqrt(1 - pre_rho^2)/(9 * pre_rho^2 + 6)/n)^(1/8)
    
    Xi <- rep(Si, each = n)
    Yi <- rep(Ti, each = n)
    Xj <- rep(Si, times = n)
    Yj <- rep(Ti, times = n)
    
    B <- rbind(2 - Xi^2 - Yi^2, mean(Si * Ti) - Xi *
                   Yi)
    pseudoB <- rbind(2 - Si^2 - Ti^2, mean(Si * Ti) -
                         Si * Ti)
    
    C1 <- ((B %*% (t(B) * matrix(rep(dnorm((Xi - Xj)/b) *  dnorm((Yi - Yj)/b), dim(B)[1]),
                                 n^2,
                                 dim(B)[1]))) - pseudoB %*% t(pseudoB) *
               dnorm(0)^2)/n/(n - 1)/b^2
    C21 <- (colSums(t(B) * matrix(rep((((Xi - Xj)/b)^2 - 1) * dnorm((Xi - Xj)/b) *
                                          dnorm((Yi - Yj)/b), dim(B)[1]),
                                  n^2,
                                  dim(B)[1])) + rowSums(dnorm(0)^2 * pseudoB))/n/(n - 1)/b^4
    C22 <- (colSums(t(B) * matrix(rep((((Yi - Yj)/b)^2 - 1) * dnorm((Xi - Xj)/b) *
                                          dnorm((Yi - Yj)/b),  dim(B)[1]),
                                  n^2,
                                  dim(B)[1])) + rowSums(dnorm(0)^2 * pseudoB))/n/(n - 1)/b^4
    
    phi40 <- sum((((Xi - Xj)/b)^4 - 6 * ((Xi - Xj)/b)^2 + 3) *
                     dnorm((Xi - Xj)/b) * dnorm((Yi - Yj)/b))/n^2/b^6
    phi04 <- sum((((Yi - Yj)/b)^4 - 6 * ((Yi - Yj)/b)^2 + 3) *
                     dnorm((Xi - Xj)/b) * dnorm((Yi - Yj)/b))/n^2/b^6
    phi22 <- sum((((Xi - Xj)/b)^2 - 1) * (((Yi - Yj)/b)^2 -  1) *
                     dnorm((Xi - Xj)/b) * dnorm((Yi - Yj)/b))/n^2/b^6
    
    if (rho.add) {
        phi31 <- sum((((Xi - Xj)/b)^3 - 3 * ((Xi -  Xj)/b)) * ((Yi - Yj)/b) *
                         dnorm((Xi - Xj)/b) * dnorm((Yi - Yj)/b))/n^2/b^6
        phi13 <- sum((((Yi - Yj)/b)^3 - 3 * ((Yi - Yj)/b)) * ((Xi - Xj)/b) *
                         dnorm((Xi - Xj)/b) *  dnorm((Yi - Yj)/b))/n^2/b^6
    } else {
        phi31 <- phi13 <- 0
    }
    
    if (rho.add) {
        C23 <- (colSums(t(B) * matrix(rep(((Yi - Yj) * (Xi - Xj)/b^2) *  dnorm((Xi - Xj)/b) *
                                              dnorm((Yi -   Yj)/b), dim(B)[1]),
                                      n^2,
                                      dim(B)[1])) + rowSums(dnorm(0)^2 * pseudoB))/n/(n - 1)/b^4
    } else {
        C23 <- matrix(rep(0, dim(B)[1]), dim(B)[1], 1)
    }
    
    if (rho.add) {
        M <- function(param) {
            rho <- 2 * pnorm(param) - 1
            C2 <- matrix(C21 + C22 + 2 * rho * C23, dim(B)[1], 1)
            C3 <- phi40 + phi04 + (4 * rho^2 + 2) *  phi22 + 4 * rho * (phi31 + phi13)
            (C3 - t(C2) %*% solve(C1) %*% C2)/(1 - rho^2)
        }
        
        rho_optimization <- optim(0,
                                  M,
                                  method = "BFGS",
                                  control = list(maxit = 20000))
        if (rho_optimization$convergence != 0) {
            stop("Check the optimization in choosing rho.")
        }
        
        rho <- 2 * pnorm(rho_optimization$par) - 1
    } else {
        rho <- 0
    }
    
    C2 <- matrix(C21 + C22 + 2 * rho * C23, dim(B)[1],
                 1)
    C3 <- phi40 + phi04 + (4 * rho^2 + 2) * phi22 +
        4 * rho * (phi31 + phi13)
    
    h <- as.numeric((1/2/pi/n/
                         (C3 - t(C2) %*% solve(C1) %*%  C2)/
                         sqrt(1 - rho^2))^(1/6))
    theta <- as.vector(-h^2/2 * solve(C1) %*% C2)
    
    c(h, rho, theta)
}

bw_tt_pcv <- function(obs, rho.add = T) {
    # This function uses the profile cross validation
    # method to select the optimal smoothing
    # parameters.
    n <- dim(obs)[1]
    Si <- qnorm(obs[, 1])
    Ti <- qnorm(obs[, 2])
    
    pre_rho <- mean(Si * Ti)
    a2 <- 19/4/pi^2
    a3 <- (48 * pre_rho^2 + 57)/32/pi^2/(pre_rho^2 -
                                             1)^3/sqrt(1 - pre_rho^2)
    a4 <- (18 * pre_rho^6 + 360 * pre_rho^4 + 576 *
               pre_rho^2 + 171)/(256 * pi^2 * (1 - pre_rho^2)^7)
    b <- (6 * a2/(sqrt(a3^2 + 3 * a2 * a4) - a3)/n)^(1/8)
    
    Xi <- rep(Si, each = n)
    Yi <- rep(Ti, each = n)
    Xj <- rep(Si, times = n)
    Yj <- rep(Ti, times = n)
    
    B <- rbind(2 - Xi^2 - Yi^2, mean(Si * Ti) - Xi *
                   Yi)
    pseudoB <- rbind(2 - Si^2 - Ti^2, mean(Si * Ti) -
                         Si * Ti)
    
    C1 <- ((B %*% (t(B) * matrix(rep(dnorm((Xi - Xj)/b) *
                                         dnorm((Yi - Yj)/b), dim(B)[1]), n^2, dim(B)[1]))) -
               pseudoB %*% t(pseudoB) * dnorm(0)^2)/n/(n -
                                                           1)/b^2
    C21 <- (colSums(t(B) * matrix(rep((((Xi - Xj)/b)^2 -  1) * dnorm((Xi - Xj)/b) *
                                          dnorm((Yi - Yj)/b),  dim(B)[1]), n^2, dim(B)[1])) +
                rowSums(dnorm(0)^2 *  pseudoB))/n/(n - 1)/b^4
    C22 <- (colSums(t(B) * matrix(rep((((Yi - Yj)/b)^2 -  1) * dnorm((Xi - Xj)/b) *
                                          dnorm((Yi - Yj)/b), dim(B)[1]), n^2, dim(B)[1])) +
                rowSums(dnorm(0)^2 * pseudoB))/n/(n - 1)/b^4
    C23 <- (colSums(t(B) * matrix(rep(((Yi - Yj) * (Xi - Xj)/b^2) *
                                          dnorm((Xi - Xj)/b) * dnorm((Yi -  Yj)/b), dim(B)[1]), n^2, dim(B)[1])) +
                rowSums(dnorm(0)^2 * pseudoB))/n/(n - 1)/b^4
    
    phi40 <- sum((((Xi - Xj)/b)^4 - 6 * ((Xi - Xj)/b)^2 +
                      3) * dnorm((Xi - Xj)/b) * dnorm((Yi - Yj)/b))/n^2/b^6
    phi04 <- sum((((Yi - Yj)/b)^4 - 6 * ((Yi - Yj)/b)^2 +
                      3) * dnorm((Xi - Xj)/b) * dnorm((Yi - Yj)/b))/n^2/b^6
    phi22 <- sum((((Xi - Xj)/b)^2 - 1) * (((Yi - Yj)/b)^2 - 1) *
                     dnorm((Xi - Xj)/b) * dnorm((Yi - Yj)/b))/n^2/b^6
    phi31 <- sum((((Xi - Xj)/b)^3 - 3 * ((Xi - Xj)/b)) * ((Yi - Yj)/b) *
                     dnorm((Xi - Xj)/b) * dnorm((Yi - Yj)/b))/n^2/b^6
    phi13 <- sum((((Yi - Yj)/b)^3 - 3 * ((Yi - Yj)/b)) * ((Xi - Xj)/b) *
                     dnorm((Xi - Xj)/b) * dnorm((Yi -  Yj)/b))/n^2/b^6
    
    wcv <- function(h, rho) {
        C2 <- matrix(C21 + C22 + 2 * rho * C23, dim(B)[1], 1)
        theta <- as.vector(solve(C1) %*% C2) * (-h^2/2)
        delta <- sqrt(h^4 * (1 - rho^2) * (4 * theta[1]^2 -  theta[2]^2) +
                          2 * h^2 * (2 * theta[1] + rho * theta[2]) + 1)
        eta <- mean(exp(-((4 * h^2 * theta[1]^2 - h^2 * theta[2]^2 + 2 * theta[1]) *
                              (Si^2 + Ti^2) +  (2 * rho * h^2 * theta[2]^2 - 8 * rho *
                                                    h^2 * theta[1]^2 + 2 * theta[2]) * Si * Ti)/
                            2/delta^2))/delta
        
        alpha1 <- -2 * h^4 * (1 - rho^2) * (4 * theta[1]^2 -  theta[2]^2) +
            2 * h^2 * ((rho^2 - 3) * theta[1] - rho * theta[2]) - 1
        alpha2 <- 4 * h^4 * (1 - rho^2) * (4 * theta[1]^2 -  theta[2]^2) +
            2 * h^2 * (4 * rho * theta[1] + (3 * rho^2 - 1) * theta[2]) + 2 * rho
        alpha3 <- -2 * h^2 * (4 * rho * theta[1] + (1 + rho^2) * theta[2]) - 2 * rho
        alpha4 <- 4 * h^2 * ((1 + rho^2) * theta[1] + rho * theta[2]) + 2
        
        part1 <- mean(exp((alpha1 * (Xi^2 + Xj^2 + Yi^2 + Yj^2) +
                               alpha2 * (Xi * Yi + Xj * Yj) +
                               alpha3 * (Xi * Yj + Xj * Yi) + alpha4 *
                               (Xi * Xj + Yi * Yj))/4/h^2/(1 - rho^2)/delta^2))/
            4/pi/eta^2/h^2/delta/sqrt(1 - rho^2)
        
        part2 <- sum(sapply(1:n,
                            function(index)
                                as.numeric(eval_tt(matrix(obs[index, ], ncol = 2),
                                                   obs[-index, ],
                                                   c(h, rho, theta))) *
                                dnorm(qnorm(obs[index, 1])) *
                                dnorm(qnorm(obs[index, 2]))
        )
        )
        part1 - 2 * part2/n
    }
    
    if (rho.add) {
        M <- function(param) {
            rho <- 2 * pnorm(param) - 1
            C2 <- matrix(C21 + C22 + 2 * rho * C23,
                         dim(B)[1], 1)
            C3 <- phi40 + phi04 + (4 * rho^2 + 2) *
                phi22 + 4 * rho * (phi31 + phi13)
            (C3 - t(C2) %*% solve(C1) %*% C2)/(1 -
                                                   rho^2)
        }
        
        rho_optimization <- optim(0,
                                  M,
                                  method = "BFGS",
                                  control = list(maxit = 20000))
        opt_rho <- 2 * pnorm(rho_optimization$par) - 1
        obj <- function(param) wcv(exp(param), opt_rho)
        optimization <- optim(log(0.2),
                              obj,
                              method = "BFGS",
                              control = list(maxit = 50000))
        opt_h <- exp(optimization$par)
    } else {
        obj <- function(param) wcv(exp(param), 0)
        optimization <- optim(log(0.2),
                              obj,
                              method = "BFGS",
                              control = list(maxit = 50000))
        opt_h <- exp(optimization$par)
        opt_rho <- 0
    }
    opt_C2 <- matrix(C21 + C22 + 2 * opt_rho * C23, dim(B)[1], 1)
    opt_theta <- as.vector(solve(C1) %*% opt_C2) * (-opt_h^2/2)
    c(opt_h, opt_rho, opt_theta, optimization$convergence)
}
