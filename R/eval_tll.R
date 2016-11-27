eval_tll <- function(uev, lfit, B) {
    uev <- as.matrix(uev)
    if(ncol(uev) == 1L) 
        uev <- matrix(uev, 1L, nrow(uev))
    d <- ncol(uev)
    zev <- qnorm(uev)
    ev  <- zev %*% solve(B)
    
    rescale <- pmax(apply(dnorm(zev), 1L, prod), 10^(- 2 * d)) * abs(det(B))
    suppressWarnings(as.numeric(predict(lfit, ev) / rescale))
}
