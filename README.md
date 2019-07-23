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

devtools::install_github("MRCIEU/TwoSampleMR") # this is a dependency not present in R CRAN, it should be installed before reverseC

devtools::install_github("SharonLutz/reverseC")
```
The install process will involve compiling source code. If you are on MacOSX, this may involve the clang compiler issuing warnings about unknown pragmas similar to the text below. Do not be alarmed if you see these. If there is actually an error, it will be present among the last several messages issued by the compiler.
```
warning: pragma diagnostic pop could not pop, no matching push [-Wunknown-pragmas]
```