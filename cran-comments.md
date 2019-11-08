## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* Suggests or Enhances not in mainstream repositories:
  INLA
  (INLA is an R package not on CRAN)
  
* Fixed ambiguous call of c++ function `pow()` under Solaris OS

* For issues from Valgrind, we carefully avoid non-initialized items. 
  The source of the uninitialized value seems to be a stack allocation in the function 
  `dlatrs_` (in `/usr/lib/libopenblasp-r0.2.19.so`), not inside `phyr` code.  