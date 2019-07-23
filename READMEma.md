

## reverseC for Mediation Analysis
These functions examines the performance of mediation analysis methods in the presence of reverse causality.

## Installation
```
you need R v3.4 or higher

You will need the proper compiling tools for your platform.
For Windows (Rtools installer):  https://cran.r-project.org/bin/windows/Rtools/
For MacOSX (clang and gfortran): https://cran.r-project.org/bin/macosx/tools/

install.packages("devtools") # devtools must be installed first

devtools::install_github("SharonLutz/reverseC")
```

## Example:
```
library(reverseC)
?reverseMAsim # For details on this function

reverseMAsim(n = 1000, pX = 0.2, gamma0 = 0, gammaX = 0.2, varM = 1, beta0 = 0, betaX = 0, 
betaM = c(0.1, 0.2, 0.3), varY = 1, nSim = 500, nSimImai = 500, SEED = 1, plot.pdf = T, 
plot.name = "reverseMAplot.pdf", alpha_level = 0.05)

 reverseMAsim(n = 1000, pX = 0.2, gamma0 = 0, gammaX = 0, varM = 1, beta0 = 0, betaX = 0.2, 
betaM = c(0.1, 0.2, 0.3), varY = 1, nSim = 500, nSimImai = 500, SEED = 1, plot.pdf = T, 
plot.name = "reverseMAplotDirect.pdf", alpha_level = 0.05)

 reverseMAsim(n = 1000, pX = 0.2, gamma0 = 0, gammaX = 0.2, varM = 1, beta0 = 0, betaX = 0.2, 
betaM = c(0.1, 0.2, 0.3), varY = 1, nSim = 500, nSimImai = 500, SEED = 1, plot.pdf = T, 
plot.name = "reverseMAplotBoth.pdf", alpha_level = 0.05)
```

## Speeding things up with optional parameters:
```
MultiProcessing And/Or Use of Threading and Eigen via C++
the reverseMAsim command accepts the following parameters:
use_multi_processing, a boolean (T, F, True, or False), which turns on the multi-processing feature
use_cpp, a boolean(T, F, True, or False), which activates the use of Rcpp RcppEigen, and threading if multiprocessing is not turned on as well.
num_jobs, an integer specifying the number of processes or threads you wish to use.

Example Using Rcpp with Eigen and 5 threads:

reverseMAsim(n = 1000, pX = 0.2, gamma0 = 0, gammaX = 0.2, varM = 1, beta0 = 0, betaX = 0, 
betaM = c(0.1, 0.2, 0.3), varY = 1, nSim = 500, nSimImai = 500, SEED = 1, plot.pdf = T, 
plot.name = "reverseMAplot.pdf", alpha_level = 0.05, use_cpp=T, num_jobs=5)


Example Using MultiProcessing and vanilla R (without using Rcpp with Eigen) with 7 subprocesses:

reverseMAsim(n = 1000, pX = 0.2, gamma0 = 0, gammaX = 0, varM = 1, beta0 = 0, betaX = 0.2, 
betaM = c(0.1, 0.2, 0.3), varY = 1, nSim = 500, nSimImai = 500, SEED = 1, plot.pdf = T, 
plot.name = "reverseMAplotDirect.pdf", alpha_level = 0.05, use_multi_processing=T, num_jobs=7)


Example Using MultiProcessing and Rcpp with Eigen with 4 subprocesses (1 thread each):

reverseMAsim(n = 1000, pX = 0.2, gamma0 = 0, gammaX = 0.2, varM = 1, beta0 = 0, betaX = 0.2, 
betaM = c(0.1, 0.2, 0.3), varY = 1, nSim = 500, nSimImai = 500, SEED = 1, plot.pdf = T, 
plot.name = "reverseMAplotBoth.pdf", alpha_level = 0.05, use_cpp=T, use_multi_processing=T, num_jobs=4)
```

## Important Caveats for the Faster Processing Strategies:
```
It is advisable that you tailor your num_jobs variable to be the # of your CPU cores - 1, at the maximum, to leave 1 core free to handle the original calling R process and any background OS processes. If you use too many procesess or threads, your system will become slow and relatively unresponsive and may lock up until the procssing completes.

The Multiprocessing strategies will use significantly more RAM, as they will store approximately the same amount of data per process that the vanilla approaches and Rcpp with Eigen and threading approach will store. If you swamp out your system RAM, paging will slow your system to a crawl, and may require a system restart to regain control.

If you force-stop an R terminal or R process that has already begun a multi-processing task, it may not be able to close all the child processes before it terminates. They will have to be stopped/killed before they stop consuming CPU and Memory and release the system resources.

```

## Output

<img src="https://github.com/SharonLutz/reverseC/blob/master/reverseMAplot.png" width="600">
