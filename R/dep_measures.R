#' Dependence measures of a `kdecop()` fit
#' 
#' Calculates several dependence measures derived from the copula density. All 
#' measures except `"blomqvist"` are computed by quasi Monte Carlo methods 
#' (see [rkdecop()].
#'
#' @param object an object of class `kdecopula`.
#' @param measures which measures to compute, see *Details*.
#' @param n_qmc the number of quasi Monte Carlo samples.
#' @param seed the seed for quasi Monte Carlo integration.
#'
#' @return A named vector of dependence measures.
#' 
#' The following measures are available:
#' \describe{
#' \item{`"kendall"`}{Kendall's \eqn{\tau}, see Nelsen (2007); computed as the 
#' sample version of a quasi Monte Carlo sample.}
#' \item{`"spearman"`}{Spearman's \eqn{\rho}, see Nelsen (2007); computed as the
#' sample version of a quasi Monte Carlo sample.}
#' \item{`"blomqvist"`}{Blomqvist's \eqn{\beta}, see Nelsen (2007); computed
#' as \eqn{4C(0.5, 0.5) - 1}.}
#' \item{`"gini"`}{Gini's \eqn{\gamma}, see Nelsen (2007); computed by quasi
#' Monte Carlo integration.}
#' \item{`"vd_waerden"`}{van der Waerden's coefficient, see Genest and Verret 
#' (2005); computed as the sample version of a quasi Monte Carlo sample.}
#' \item{`"minfo"`}{mutual information, see Joe (1989); computed by quasi Monte
#' Carlo integration.}
#' \item{`"linfoot"`}{Linfoot's correlation coefficient, see Joe (1989); 
#' computed by quasi Monte Carlo integration.}
#' }
#' 
#' @references 
#' Nelsen, R. (2007). An introduction to copulas. Springer Science 
#' & Business Media, 2007.  
#' 
#' Genest, C., and Verret, F. (2005). Locally most powerful rank tests of 
#' independence for copula models. Journal of Nonparametric Statistics, 17(5)  
#' 
#' Joe, H. (1989). Relative Entropy Measures of Multivariate Dependence. 
#' Journal of the American Statistical Association, 84(405)  
#'
#' @examples
#' ## load data and transform with empirical cdf
#' data(wdbc)
#' udat <- apply(wdbc[, -1], 2, function(x) rank(x) / (length(x) + 1))
#' 
#' ## estimate copula density and calculate dependence measures
#' fit <- kdecop(udat[, 5:6])
#' dep_measures(fit)
#' 
#' @export
dep_measures <- function(object, measures = "all", n_qmc = 10^3, seed = 5) {
    if (any(measures == "all")) measures <- all_measures
    
    # quasi Monte Carlo samples
    stopifnot(n_qmc >= 10)
    set.seed(seed)
    u_qmc <- rkdecop(n_qmc, object)
    
    # calculate measures
    sapply(measures, calculate_dep_measure, u = u_qmc, object = object)
}

all_measures <- c("kendall", "spearman", "blomqvist", "gini", "vd_waerden",
                  "minfo", "linfoot")

calculate_dep_measure <- function(measure, u, object) {
    result <- switch(
        measure,
        "kendall"    = stats::cor(u, method = "kendall")[1, 2],
        "spearman"   = stats::cor(u, method = "spearman")[1, 2],
        "blomqvist"  = 4 * pkdecop(c(0.5, 0.5), object) - 1,
        "gini"       = 2 * mean(abs(u[, 1] + u[, 2] - 1) - abs(u[, 1] - u[, 2])),
        "vd_waerden" = stats::cor(qnorm(u))[1, 2],
        "minfo"      = mean(log(dkdecop(u, object))),
        "linfoot"    = sqrt(1 - exp(-2 * mean(log(dkdecop(u, object)))))
    )
}
