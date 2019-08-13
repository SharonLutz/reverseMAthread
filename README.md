## reverseC
This package examines the performance of Mendelian Randomization (MR) and Mediation Analysis methods in the presence of reverse causality. 

## reverseC for Mendelian Randomization
A full a description of the R functions that examine the role of the exposure on the outcome when the outcome is correctly and incorrectly specified for 3 common MR approaches is given [here](READMEmr.md).

## reverseC for Mediation Analysis
A full a description of the R functions that examine the role of the exposure on the outcome through the mediator when the outcome is correctly and incorrectly specified for mediation analysis is given [here](READMEma.md).

## Installation
Requirements:
* R v3.4 or higher
* You will need the proper compiling tools for your platform.
  * For Windows (Rtools installer): https://cran.r-project.org/bin/windows/Rtools/
  * For MacOSX (clang and gfortran): https://cran.r-project.org/bin/macosx/tools/

```
install.packages("devtools") # devtools must be installed first

# these will fail to install when already loaded, and install_github will sometimes 
# load these as part of its activity, and will then try to install them if they need 
# an update for one of the package dependencies
install.packages(c("Rcpp","RcppEigen", "curl"), quiet=T) 

#this package does not install automatically, but is needed by TwoSampleMR
install.packages("psych")

devtools::install_github("MRCIEU/TwoSampleMR") # this is a dependency not present in R CRAN, it should be installed before reverseC

devtools::install_github("SharonLutz/reverseC",quiet=T)
```
The install process will involve compiling source code. If you are on MacOSX, this may involve the clang compiler issuing warnings about unknown pragmas similar to the text below. Do not be alarmed if you see these. If there is actually an error, it will be present among the last several messages issued by the compiler.
```
warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
```
## Warning: Do not try to access package internals directly or do so at your own risk!
If you try to run methods/functions that are not exported and intended for end users, and feed these functions environments, parameters, or values that are not correctly formed, it could result in an uncaught or uncatchable C++ exception or segmentation fault. If this occurs, it will kill your R session/terminal and if you were working within RStudio it will probably crash too.