### bandwidth selection
bw_select <- function(udata, method) {
    
    switch(method,
           "MR"   = bw_mr(udata),
           "beta" = bw_beta(udata),
           "T"    = nrow(udata)^(-1/(ncol(udata) + 4)) * t(chol(cov(qnorm(udata)))),
           "TLL1" = bw_tll(qnorm(udata), deg = 1),
           "TLL2" = bw_tll(qnorm(udata), deg = 2),
           "TLL1c" = bw_tllc(qnorm(udata), deg = 1),
           "TLL2c" = bw_tllc(qnorm(udata), deg = 2),
           "TTPI" = bw_tt_plugin(udata),
           "TTCV" = bw_tt_pcv(udata))
}

# precalculated integral values
tau.sq <- seq(-0.9, 0.9, l = 50)
beta.sq <- c(5.8e+06, 2.4e+06, 1.7e+06, 4.5e+05, 2.7e+05, 1.7e+05, 1e+05, 7e+04,
             6.3e+04, 4.9e+04, 3.1e+04, 2.1e+04, 1.5e+04, 1.1e+04, 8e+03, 5.8e+03,
             4.2e+03, 3.1e+03, 2.3e+03, 1.8e+03, 1.3e+03, 1e+03, 7.7e+02, 5.8e+02,
             4.4e+02, 3.4e+02, 2.6e+02, 2e+02, 1.5e+02, 1.1e+02,  86,  64,  48,
             35,  26,  19,  13, 9.3, 6.4, 4.3, 2.7, 1.7,   1, 0.55, 0.28, 0.12,
             0.044, 0.011, 0.0015, 1.8e-05, 1.8e-05, 0.0015, 0.011, 0.044, 0.12, 
             0.28, 0.55,   1, 1.7, 2.7, 4.3, 6.4, 9.3,  13,  19,  26,  35,  48, 
             64,  86, 1.1e+02, 1.5e+02, 2e+02, 2.6e+02, 3.4e+02, 4.4e+02, 5.8e+02, 
             7.7e+02, 1e+03, 1.3e+03, 1.8e+03, 2.3e+03, 3.1e+03, 4.2e+03, 5.8e+03, 
             8e+03, 1.1e+04, 1.5e+04, 2.1e+04, 3.1e+04, 4.9e+04, 6.3e+04, 7e+04, 
             1e+05, 1.7e+05, 2.7e+05, 4.5e+05, 1.7e+06, 2.4e+06, 5.8e+06)
xi.sq <- c(4.3e+04, 2e+04, 1.3e+04, 4.9e+03, 2.8e+03, 1.7e+03, 1.5e+03, 1.1e+03,
           7e+02, 5.1e+02, 3.7e+02, 2.9e+02, 2.2e+02, 1.7e+02, 1.3e+02, 1.1e+02,
           83,  68,  54,  45,  38,  32,  27,  23,  20,  17,  14,  12,  10,   9,
           7.7, 6.5, 5.6, 4.7,   4, 3.3, 2.8, 2.3, 1.9, 1.5, 1.2, 0.95, 0.72, 
           0.53, 0.38, 0.25, 0.15, 0.075, 0.027, 0.003, 0.003, 0.027, 0.075, 
           0.15, 0.25, 0.38, 0.53, 0.72, 0.95, 1.2, 1.5, 1.9, 2.3, 2.8, 3.3,  
           4, 4.7, 5.6, 6.5, 7.7,   9,  10,  12,  14,  17,  20,  23,  27,  32, 
           38,  45,  54,  68,  83, 1.1e+02, 1.3e+02, 1.7e+02, 2.2e+02, 2.9e+02,
           3.7e+02, 5.1e+02, 7e+02, 1.1e+03, 1.5e+03, 1.7e+03, 2.8e+03, 4.9e+03,
           1.3e+04, 2e+04, 4.3e+04)
zeta.sq <- c( 18,  21,  21,  21,  20,  20,  19,  19,  18,  18,  17,  17,  15, 
              16,  15,  15,  15,  14,  14,  14,  14,  13,  13,  13,  12,  12,  12,  12,  12, 
              11,  11,  11,  11,  11,  11,  10,  10,  10,  10, 9.9, 9.8, 9.7, 9.6, 9.5, 9.4, 
              9.4, 9.3, 9.3, 9.3, 9.3, 9.3, 9.3, 9.3, 9.3, 9.4, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9,  
              10,  10,  10,  10,  11,  11,  11,  11,  11,  11,  12,  12,  12,  12,  12,  13, 
              13,  13,  14,  14,  14,  14,  15,  15,  15,  16,  15,  17,  17,  18,  18,  19, 
              19,  20,  20,  21,  21,  21,  18)

bw_mr <- function(udata) {
    n <- nrow(udata)
    
    ## constants for kernel
    sigma_K <- sqrt(1/5)
    d_K     <- 3/5
    
    ## parameter for frank copula by inversion of Kendall's tau
    tau <- cor(udata, method="kendall")[1L, 2L]
    
    ## integrals
    tau.ind <- which.min(tau.sq - tau)
    beta  <- beta.sq[tau.ind]
    gamma <- 1
    
    ## result
    res <- (2*d_K^2/sigma_K^4*gamma/beta)^(1/6) * n^(-1/6)
    if (res > 1) 1 else res
}



bw_beta <- function(udata) {
    n  <- nrow(udata)
    
    tau <- cor(udata, method="kendall")[1L, 2L]
    
    ## integrals
    tau.ind <- which.min(tau.sq - tau)
    xi   <- xi.sq[tau.ind]
    zeta <- zeta.sq[tau.ind]
    
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
    pca <- princomp(zdata)
    B <- unclass(pca$loadings)
    B <- B * sign(diag(B))
    qrs <- unclass(pca$scores)
    
    ## find optimal alpha in for each principal component
    alphsq <- seq(nrow(zdata)^(-1/5), 1, l = 50)
    opt <- function(i) {
        val <- lscvplot(~qrs[, i], 
                        alpha  = alphsq, 
                        deg    = deg, 
                        kern   = "gauss", 
                        maxk   = 512)$value
        mean(alphsq[which.min(val)])
    }
    alpha.vec <- sapply(1:d, opt)

    ## adjustments for multivariate estimation and transformation
    kappa <- alpha.vec[1]/alpha.vec
    dimnames(B) <- NULL
    if (deg == 1) {
        alpha <- n^(1/5 - d/(4 + d)) * alpha.vec[1]
    } else {
        alpha <- n^(1/9 - d/(8 + d)) * alpha.vec[1]
    }
    
    ## return results
    list(B = B, alpha = alpha, kappa = kappa)
}

bw_tllc <- function(zdata, deg) {
    n <- nrow(zdata)
    d <- ncol(zdata)
    # transform to uncorrelated data
    H <- cov(zdata) * nrow(zdata)^(-2/ (2 * 2 * deg + d))
    
    solve(chol(H)) 
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
