#' Calculates a Bernstein polynomial
#'
#' @param x vector of evaluation points
#' @param k index of the Bernstein polynomial
#' @param m order of the Bernstein polynomial
#' @noRd
bern_poly <- function(x, k, m) {
    (choose(m, k) * x^k * (1 - x)^(m - k)) * (m + 1)
}

#' Calculates basis coefficients fo the Bernstein estimator
#'
#' The method of Weiss and Scheffer (2016) is used to adjust the coefficients
#' to uniform margins.
#'
#' @param u data
#' @param m order of the Bernstein polynomial
#' @noRd
#' @importFrom quadprog solve.QP
bern_coefs <- function(u, m) {
    # initialize coefficients with empirical frequencies
    ucut <- lapply(1:2, function(i) cut(u[, i], 0:(m + 1) / (m + 1)))
    cf0 <- table(ucut[[1]], ucut[[2]]) / nrow(u)
    as.matrix(cf0)
    # # set up quadratic programming problem
    # m <- m + 1
    # D <- diag(m^2)
    # d <- as.vector(cf0)
    # A1 <- A2 <- matrix(0, m, m^2)
    # for (i in 1:m) {
    #     A1[i, (i - 1) * m + 1:m] <- 1
    #     A2[i, m * (1:m) - m + i] <- 1
    # }
    # A3 <- diag(m^2)
    # A <- t(rbind(A1, A2, A3))
    # b <- c(rep(1 / m, 2 * m), rep(0, m^2))
    # 
    # # solve
    # sol <- tryCatch(solve.QP(D, d, A, b, meq = 2 * m)$solution,
    #                 error = function(e) d)
    # sol[sol < 0] <- 0
    
    # return as matrix
    # matrix(sol, m, m)
}

#' Fit a Bernstein copula to the data
#'
#' @param u data
#' @param m order of the Bernstein copula
#'
#' @return berncop object
#' @noRd
berncop <- function(u, m = 10) {
    out <- list(coefs = bern_coefs(u, m),  m = m)
    class(out) <- "berncop"
    out
}


#' Evaluate the density of a berncop object
#'
#' @param unew data
#' @param object berncop object
#'
#' @return Bernstein copula density evaluated at uev.
#' @noRd
dberncop <- function(uev, object) {
    if (NCOL(uev) == 1)
        uev <- matrix(uev, ncol = 2)
    list2env(object)
    summand <- function(i, j) {
        P1 <- bern_poly(uev[, 1], i, object$m)
        P2 <- bern_poly(uev[, 2], j, object$m)
        object$coefs[i + 1, j + 1] * P1 * P2
    }
    
    est <- 0
    for (i in 0:object$m) {
        phi_i <- sapply(0:object$m, function(j) summand(i, j))
        if (NCOL(phi_i) == 1)
            phi_i <- t(phi_i)
        est <- est + rowSums(phi_i)
    }
    est
}

#' Simulates from 
#'
#' @param unew data
#' @param object berncop object
#'
#' @return Bernstein copula density evaluated at uev.
#' @importFrom stats rbeta
#' @noRd
rberncop <- function(n, object) {
    m <- object$m
    idx <- sample(seq_len((m + 1)^2), n, replace = TRUE, prob = c(object$coefs))
    ij <- expand.grid(0:m, 0:m)[idx, ]
    cbind(rbeta(n, ij[, 1] + 1, m + 1 - ij[, 1]), 
          rbeta(n, ij[, 2] + 1, m + 1 - ij[, 2]))
}
