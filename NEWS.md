kdecopula (under development
-------------------------------

NEW FEATURES

   * improved bandwidth selection.
   
   * `predict()`, `fitted()`, and `simulate()` generics.
   
   * method `TLL2nn` is the default method again.


kdecopula 0.8.1
-------------------------------

NEW FEATURES

   * stabilized bandwidth matrix transformation.
   
   * faster builds (removed unnecessarily included headers).
   
   * updated documentation.


kdecopula 0.8.0
-------------------------------

NEW FEATURES

  * Improved bandwidth selection heuristic for methods `TLL1/TLL2`.
  
  * The Bernstein copula density estimator is documented and comes with an
    automated smoothing parameter selection routine.
    
  * All bandwidth selection methods now have a manual page (see `?kdecop` for
    links).
    
  * Updated vignette.

BUG FIXES

  * fixed inconsistencies in bandwidth specification. The interpretation of the
    bandwidth matrix is now the same for `TLL` and `T` methods.


kdecopula 0.7.0
-------------------------------

DEPENDS

  * Added `quadprog` to imports.

NEW FEATURES
  
  * Methods `TLL1` and `TLL2` now use fixed bandwidths; nearest neighbor methods
    are implemented as `TLL1nn` and `TLL2n`.
  
  * Improved bandwidth selection for methods `TLL1nn` and `TLL2nn`.
  
  * Improved memory allocation.
  
  * Hidden (and undocumented) method `bern` for the Bernstein copula density.
  
BUG FIXES

  * Bandwidth selection for methods `TTxx` are now safe even when correlation is
    high and sample size small.


kdecopula 0.6.0
-------------------------------

DEPENDS

  * Package no longer imports packages `cubature` and `VineCopula`.

BUG FIXES

  * improved algorithms for numerical inversion
  
NEW FEATURES

  * new function `hkdecop` to evaluate conditional distributions 
    ("h-functions").
    
  * extended plotting functionality (exponential margins).
  
  * faster calculations of bandwidths for methods "MR" and "beta" due to
    precalculated integral values.
    
OTHER

  * updated vignette.
  
  
kdecopula 0.5.0
-------------------------------

BUG FIXES

  * fix scaling for method = "MR" when not renormalized.
  
NEW FEATURES

  * package now has a vignette.



kdecopula 0.4.1
-------------------------------

BUG FIXES

  * rkdecop now properly working with quasi-random numbers.
 
  * improved efficiency of C++ code. 


kdecopula 0.4.0
-------------------------------

NEW FEATURES

  * New estimation method: tapered transformation estimator 
    (methods "TTPI" and "TTCV"). 
    Thanks to Kuangyu Wen for sharing code!

kdecopula 0.3.0
-------------------------------

NEW FEATURES

  * New methods for class `kdecopula`: `logLik`, `AIC`, `BIC`
  
  * Full summary (logLik, AIC, BIC, effp) even when `info = FALSE`.
  

BUG FIXES

  * Fixed typo in install instructions.
  
  * Fixed bug in summary (`method = "T"`).

  * Improve formatting in summary.


kdecopula 0.2.0
-------------------------------

NEW FEATURES

  * New methods for class `kdecopula`: `print`, `summary`, `contour`.
  
  * `kdecop`: Default for `info` argument is now `TRUE`.
  
  * Updated documentation.

BUG FIXES

  * `kdecop`: fix error when requesting info about fit. 


kdecopula 0.1.0
-------------------------------

DEPENDS

  * Package now imports the qrng pacakge.
  
NEW FEATURES

  * `rkdecop`: add quasi-random number option based on package qrng.


kdecopula 0.0.4
-------------------------------

BUG FIXES

   * `dkdecop`: set cutoff for `stable` option in  to 50
   
   * Fix memory read errors.


kdecopula 0.0.3
-------------------------------

BUG FIXES

  * Dixed some memory issues in C++ code.
  
  * Adjust imports for R-3.3.0.


kdecopula 0.0.2
-------------------------------

BUG FIXES

  * Fixed issues with non-portable C++ code (problems occured on solaris).
