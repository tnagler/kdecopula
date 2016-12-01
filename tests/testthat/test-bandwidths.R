library(VineCopula) 
context("Bandwidths") 
test_that("selection works for method 'bern'", {
    set.seed(5) 
    u1 <- cbind(runif(30), runif(30)) # low dependence 
    u2 <- cbind(u1[, 1], 1 - u1[, 2]) # low dependence (opposite sign) 
    u3 <- BiCopSim(30, 3, 10)         # strong positive dependence 
    u4 <- cbind(u3[, 1], 1 - u3[, 2]) # strong negative dependence 
    
    expect_gt1 <- function(bw) expect_gt(bw, 1) 
    expect_gt1(kdecopula:::bw_bern(u1)) 
    expect_gt1(kdecopula:::bw_bern(u2)) 
    expect_gt1(kdecopula:::bw_bern(u3)) 
    expect_gt1(kdecopula:::bw_bern(u4)) 
    
    expect_integer <- function(bw) expect_identical(bw, round(bw)) 
    expect_integer(kdecopula:::bw_bern(u1)) 
    expect_integer(kdecopula:::bw_bern(u2)) 
    expect_integer(kdecopula:::bw_bern(u3)) 
    expect_integer(kdecopula:::bw_bern(u4)) 
})
