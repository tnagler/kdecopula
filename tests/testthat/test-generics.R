context("S3 generics")

data(wdbc)
u <- VineCopula::pobs(wdbc[1:20, 6:7])
est <- kdecop(u, method = "TLL1nn")

test_that("print and summary work without error", {
    expect_output(print(est))
    expect_output(summary(est))
})

test_that("predict calls correct (d/p/h)kdecop function", {
    expect_equal(predict(est, u), dkdecop(u, est))
    expect_equal(predict(est, u, "cdf"), pkdecop(u, est))
    expect_equal(predict(est, u, "hfunc1"), hkdecop(u, est, 1))
    expect_equal(predict(est, u, "hfunc2"), hkdecop(u, est, 2))
    expect_equal(predict(est, u, "hinv1"), hkdecop(u, est, 1, inverse = TRUE))
    expect_equal(predict(est, u, "hinv2"), hkdecop(u, est, 2, inverse = TRUE))
})

test_that("fitted equals predict on original data", {
    expect_equal(predict(est, u), fitted(est))
    expect_equal(predict(est, u, "cdf"), fitted(est, "cdf"))
    expect_equal(predict(est, u, "hfunc1"), fitted(est, "hfunc1"))
    expect_equal(predict(est, u, "hfunc2"), fitted(est, "hfunc2"))
    expect_equal(predict(est, u, "hinv1"), fitted(est, "hinv1"))
    expect_equal(predict(est, u, "hinv2"), fitted(est, "hinv2"))
})

test_that("simulate equals rkdecop output", {
    set.seed(1)
    sim <- rkdecop(500, est)
    set.seed(1)
    simq <- rkdecop(500, est, quasi = TRUE)
    expect_equal(simulate(est, 500, seed = 1), sim)
    expect_equal(simulate(est, 500, seed = 1, quasi = TRUE), simq)
})
