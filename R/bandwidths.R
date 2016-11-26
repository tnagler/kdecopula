bw_select <- function(udata, method) {
    switch(method,
           "MR"     = bw_mr(udata),
           "beta"   = bw_beta(udata),
           "T"      = bw_t(udata),
           "TLL1nn" = bw_tllnn(udata, deg = 1),
           "TLL2nn" = bw_tllnn(udata, deg = 2),
           "TLL1"   = bw_tll(udata, deg = 1),
           "TLL2"   = bw_tll(udata, deg = 2),
           "TTPI"   = bw_tt_plugin(udata),
           "TTCV"   = bw_tt_pcv(udata),
           "bern"   = bw_bern(udata))
}


bw_bern <- function(udata) {
    n <- nrow(udata)
    rho <- cor(udata)[1, 2]
    # Rose (2015)
    n^(1/3) * exp(rho^(1/n)) * (rho + 0.1)
}

## mirror reflection ------------------------
bw_mr <- function(udata) {
    # optimal bandwidth depends only on |tau| 
    tau <- abs(cor(udata, method = "kendall")[1L, 2L])
    res <- precalc_bw_mr(tau) * nrow(udata)^(-1/6)
    if (res > 1) 1 else res
}



precalc_bw_mr <- function(tau) {
    # ## constants for kernel
    # sigma_K <- sqrt(1/5)
    # d_K     <- 3/5
    # 
    # ## parameter for frank copula by inversion of Kendall's tau
    # family <- 5
    # if (abs(tau) < 1e-2) {
    #     family <- 0
    #     par <- 0
    # } else {
    #     par <- BiCopTau2Par(family, tau = tau)
    # }
    # 
    # ## short handles for copula density and derivatives
    # c_uu <- function(u,v)
    #     BiCopDeriv2(u, v, family, par, deriv = "u1")
    # c_vv <- function(u,v)
    #     BiCopDeriv2(u, v, family, par, deriv = "u2")
    # 
    # ## integrals
    # require(cubature)
    # bet <- function(w) (c_uu(w[1L], w[2L]) + c_vv(w[1L], w[2L]))^2
    # beta  <- adaptIntegrate(bet,
    #                         lowerLimit = c(0, 0),
    #                         upperLimit = c(1, 1),
    #                         tol = 5e-3)$integral
    # gamma <- 1
    # 
    # ## result
    # (2*d_K^2/sigma_K^4*gamma/beta)^(1/6)
    
    ## the above calculations were done on a grid for Kendall's tau.
    v <- c(Inf, Inf, Inf, Inf, Inf, 9.38981068252189, 7.36254398988451,
           5.9934473841143, 5.01438723715381, 4.28489907652576, 3.72151713690116,
           3.27524808116924, 2.91396051165559, 2.61608148151527, 2.36661799041074,
           2.15485773318137, 1.97295862876476, 1.81504529817913, 1.67664440477963,
           1.55428808465057, 1.44524733239357, 1.34727336943063, 1.2588505869573,
           1.17831883935096, 1.10458837918203, 1.03668522552123, 0.973750351604238,
           0.915194628199262, 0.860406562943529, 0.808890181007853, 0.760230558570832,
           0.714051901933609, 0.670076097370562, 0.628067045708522, 0.587716317232067,
           0.548818087457592, 0.51120935335435, 0.474708408418336, 0.439129770121076,
           0.404297062441312, 0.370046576741511, 0.336198423712505, 0.302561619578939,
           0.26895782427755, 0.23518179071423, 0.201000830641606, 0.16613394170313,
           0.13018932667288, 0.0925244955920096, 0.037769121689797)
    ## we choose the value that correponds to the value of Kendall's tau that
    ## is closest to the empirical one.
    tausq <- seq(0, 0.98, l = 50)^2
    v[which.min(abs(tau - tausq))]
}


## beta kernels -----------------------------
bw_beta <- function(udata) {
    # optimal bandwidth depends only on |tau| 
    tau <- abs(cor(udata, method="kendall")[1L, 2L])
    precalc_bw_beta(tau) * nrow(udata)^(-1/3)
}

precalc_bw_beta <- function(tau) {
    # ## parameter for frank copula by inversion of Kendall's tau
    # family <- 5
    # if (abs(tau) < 1e-2) {
    #     family <- 0
    #     par <- 0
    # } else {
    #     par <- BiCopTau2Par(family, tau = tau)
    # }
    # 
    # ## short handles for copula density and derivatives
    # cd   <- function(u,v)
    #     BiCopPDF(u, v, family, par)
    # c_u  <- function(u,v)
    #     BiCopDeriv(u, v, family, par, deriv = "u1")
    # c_v  <- function(u,v)
    #     BiCopDeriv(u, v, family, par, deriv = "u2")
    # c_uu <- function(u,v)
    #     BiCopDeriv2(u, v, family, par, deriv = "u1")
    # c_vv <- function(u,v)
    #     BiCopDeriv2(u, v, family, par, deriv = "u2")
    # 
    # ## short handles for integrands
    # x <- function(w) {
    #     u <- w[1L]
    #     v <- w[2L]
    #     ((1-2*u)*c_u(u,v) + (1-2*v)*c_v(u,v) +
    #             1/2 * (u*(1-u)*c_uu(u,v) + v*(1-v)*c_vv(u,v)))^2
    # }
    # zet <- function(w) {
    #     u <- w[1L]
    #     v <- w[2L]
    #     cd(u,v) / sqrt(u*(1-u)*v*(1-v))
    # }
    # 
    # ## integrations
    # require(cubature)
    # xi   <- adaptIntegrate(x,
    #                        lowerLimit = c(0, 0),
    #                        upperLimit = c(1, 1),
    #                        tol = 5e-3,
    #                        maxEval = 10^3)$integral
    # zeta <- adaptIntegrate(zet,
    #                        lowerLimit = c(0, 0),
    #                        upperLimit = c(1, 1),
    #                        tol = 5e-3,
    #                        maxEval = 10^3)$integral
    # 
    # ## result
    # (zeta/(8*pi*xi))^(1/3)
    
    ## the above calculations were done on a grid for Kendall's tau.
    v <- c(Inf, Inf, Inf, Inf, Inf, 4.69204254128018, 3.67817865896611,
           2.99306705888, 2.50600489245533, 2.14070935782172, 1.85882498366308,
           1.63542027413745, 1.45447243427281, 1.30516553259084, 1.18006222840396,
           1.07369616223661, 0.982279653155062, 0.902786075906722, 0.832989698639277,
           0.77116216500065, 0.71594673568386, 0.666226838391616, 0.621283273241791,
           0.580285687536092, 0.542978050804531, 0.508272377148461, 0.476077286848213,
           0.445778127751307, 0.417823893972421, 0.391328969202865, 0.366246870917148,
           0.342361362868933, 0.31955865730248, 0.297652278179634, 0.276344663543029,
           0.255870283604903, 0.235590782152833, 0.215830541690526, 0.196004036954186,
           0.176279793709377, 0.156562165388165, 0.13363202871254, 0.116793384758542,
           0.0963960312149928, 0.0805978666127558, 0.0612707894515751, 0.0402224124093885,
           0.0294236043019451, 0.0162264027190233, 0.00459925701909334)
    ## we choose the value that correponds to the value of Kendall's tau that
    ## is closest to the empirical one.
    tausq <- seq(0, 0.98, l = 50)^2
    v[which.min(abs(tau - tausq))]
    
}

## transformation kernel -----------------------------
bw_t <- function(udata) {
    ## normal rederence rule
    n <- nrow(udata)
    n^(-1 / 6) * t(chol(cov(qnorm(udata))))
}

## transformation local likelihood -------------------
bw_tll <- function(udata, deg) {
    n <- nrow(udata)
    # transform to uncorrelated data
    sqrt(5)^deg * n^(-1 / (4 * deg + 2)) * t(chol(cov(qnorm(udata))))
}

bw_tllnn <- function(udata, deg) {
    zdata <- qnorm(udata)
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

## tapered transformation estimator ------------------
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
    # standardization is asymptotically negligble, but prevents invalid choices
    # of the bandwidth parameters when pre_rho^2 > 1
    Si <- Si / sd(Si)
    Ti <- Ti / sd(Ti)

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
                                  control = list(maxit = 2000))
        if (rho_optimization$convergence != 0) {
            stop("Check the optimization in choosing rho.")
        }
        
        rho <- 2 * pnorm(rho_optimization$par) - 1
    } else {
        rho <- 0
    }

    C2 <- matrix(C21 + C22 + 2 * rho * C23, dim(B)[1], 1)
    C3 <- phi40 + phi04 + (4 * rho^2 + 2) * phi22 + 4 * rho * (phi31 + phi13)
    h <- as.numeric((1/2/pi/n/(C3 - t(C2) %*% solve(C1) %*%  C2)/
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
    # standardization is asymptotically negligble, but prevents invalid choices
    # of the bandwidth parameters when pre_rho^2 > 1
    Si <- Si / sd(Si)
    Ti <- Ti / sd(Ti)
    
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
            (C3 - t(C2) %*% solve(C1) %*% C2)/(1 - rho^2)
        }
        
        rho_optimization <- optim(0,
                                  M,
                                  method = "BFGS",
                                  control = list(maxit = 2000))
        opt_rho <- 2 * pnorm(rho_optimization$par) - 1
        obj <- function(param) wcv(exp(param), opt_rho)
        optimization <- optim(log(0.2),
                              obj,
                              method = "BFGS",
                              control = list(maxit = 2000))
        opt_h <- exp(optimization$par)
    } else {
        obj <- function(param) wcv(exp(param), 0)
        optimization <- optim(log(0.2),
                              obj,
                              method = "BFGS",
                              control = list(maxit = 2000))
        opt_h <- exp(optimization$par)
        opt_rho <- 0
    }
    opt_C2 <- matrix(C21 + C22 + 2 * opt_rho * C23, dim(B)[1], 1)
    opt_theta <- as.vector(solve(C1) %*% opt_C2) * (-opt_h^2/2)
    opt_theta[1] <- max(opt_theta[1], 0)
    c(opt_h, opt_rho, opt_theta, optimization$convergence)
}
