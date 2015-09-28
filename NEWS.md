kdecopula 0.3.0
-------------------------------

NEW FEATURES

  * New methods for class `kdecopula`: `logLik`, `AIC`, `BIC`
  
  * Full summary (logLik, AIC, BIC, effp) even when `info = FALSE`.
  

BUG FIXES

  * Fixed type in install instructions.
  
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
