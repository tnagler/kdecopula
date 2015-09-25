<!-- README.md is generated from README.Rmd. Please edit that file -->
kdecopula
=========

> Kernel smoothing for bivariate copula densities

[![Build status Linux](https://travis-ci.org/tnagler/kdecopula.svg?branch=master)](https://travis-ci.org/tnagler/kdecopula) [![Windows Build status](http://ci.appveyor.com/api/projects/status/github/tnagler/kdecopula?svg=true)](https://ci.appveyor.com/project/tnagler/kdecopula) [![CRAN version](http://www.r-pkg.org/badges/version/kdecopula)](https://cran.r-project.org/web/packages/kdecopula/index.html) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/kdecopula)](https://cran.r-project.org/web/packages/kdecopula/index.html)

This package provides fast implementations of kernel estimators for the copula density. Due to its several plotting options it is particularly useful for the exploratory analysis of dependece structures. It can be further used for flexible nonparametric estimation of copula densities and resampling.

You can install:

-   the stable release on CRAN:

    ``` r
    install.pacakges("VineCopula")
    ```

-   the latest development version:

    ``` r
    devtools::install_github("tnagler/VineCopula")
    ```

------------------------------------------------------------------------

Functions
---------

The package provides the following functions:

-   `kdecop`: Kernel estimation of a copula density. By default, estimation method and bandwidth are selected automatically. Returns an object of class `kdecopula`.

-   `dkdecop`: Evaluates the density of a `kdecopula` object.

-   `pkdecop`: Evaluates the distribution function of a `kdecopula` object.

-   `rkdecop`: Simulates synthetic data from a `kdecopula` object.

-   Methods for calss `kdecopula`:
    -   `plot`, `contour`: Surface and contour plots of the
    -   `print`, `summary`: Displays further information about the density estimate.

Look up the package documentation for more details on arguments and options,

------------------------------------------------------------------------

kdecopula in action
-------------------

In this document, we demonstrate the main capabilities of the `kdecopula` package. All user-level functions will be introduced and demonstrated on simulated data.

Let's simulate data from a Clayton copula via `BiCopSim` from the `VineCopula` package.

``` r
library(VineCopula)
clay3 <- BiCop(family = 3, par = 3)
u <- BiCopSim(500, clay3)
plot(u)
```

![](inst/README-unnamed-chunk-3-1.png)

#### Estimation of bivariate copula densities

We start by estimating the copula density with the `kdecop` function. There is a number of options for the smoothing parameterization, estimation method and evaluation grid, but it is only required to provide a data-matrix.

``` r
kde.fit <- kdecop(u)
summary(kde.fit)
#> Kernel copula density estimate (class 'kdecopula') 
#> -------------------------------------------------- 
#> Observations: 500 
#> Method:       TLL2 
#> Bandwidth:    alpha = 0.5558749
#>               B = matrix(c(0.69, 0.72, -0.5, 0.48), 2, 2)
#> -------------------------------------------------- 
#> logLik: 294.84    AIC: -563.57    cAIC: -562.82    BIC: -508.56 
#> Effective number of parameters: 13.05
```

The output of the function `kdecop` is an object of class `kdecopula` that contains all information collected during the estimation process and summary statistics such as *AIC* or the *effective number of parameters*.

#### Working with a `kdecopula` object

The density and *cdf* can be computed easily:

``` r
dkdecop(c(0.1, 0.2), kde.fit)
#> [1] 2.384852
pkdecop(c(0.1, 0.2), kde.fit)
#> [1] 0.09245737
```

Furthermore, we can simulate synthetic data from the estimated density:

``` r
unew <- rkdecop(500, kde.fit)
plot(unew)
```

![](inst/README-unnamed-chunk-6-1.png)

#### Plotting bivariate copula densities

The most interesting part for most people is probably to make exploratory plots. The class `kdecopula` has its own generic for plotting. In general, there are two possible types of plots: *contour* and *surface* (or perspective) plots. Additionally, the `margins` argument allows to choose between a plot of the original copula density and a meta-copula density with standard normal margins (default for `type = surface`).

``` r
plot(kde.fit)
```

![](inst/README-unnamed-chunk-7-1.png)

``` r
contour(kde.fit)
```

![](inst/README-unnamed-chunk-8-1.png)

``` r
contour(kde.fit, margins = "unif")
```

![](inst/README-unnamed-chunk-9-1.png)
