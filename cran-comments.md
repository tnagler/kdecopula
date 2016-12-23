## Resubmission

Problem : 
Package has a VignetteBuilder field but no prebuilt vignette index. 
* checking files in 'vignettes' ... WARNING
Files in the 'vignettes' directory newer than all files in 'inst/doc':
  'kdecopula.pdf.asis' 
  
Solution: 
Built pdf from scratch and include in `vignettes` and `inst/doc`. An index 
should be created by 'kdecopula.pdf.asis'.

## Test environments
* ubuntu 12.04 (devel, release) 
* win-builder (release)

## R CMD check results
There were no ERRORs or WARNINGs. 

## Downstream dependencies
Checked VineCopula: 0 errors | 0 warnings | 0 notes
