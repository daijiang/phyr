## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

  Suggests or Enhances not in mainstream repositories:
  INLA
  (INLA is an R package not on CRAN)
  
* We have unwrapped examples that took ~5 seconds out of \dontrun.
  Examples took >> 5 seconds to finish (phylogeneic generalized linear 
  mixed models, bootstrapping) are still wrapped within \dontrun{}.
  
* We have removed 'with R' out of the title. We have also corrected all
  T/F to TRUE/FALSE.