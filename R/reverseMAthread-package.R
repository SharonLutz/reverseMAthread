## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib reverseMAthread, .registration = TRUE
#' @import RcppEigen

## usethis namespace: end
#' @title reverseMAthread-package
#' @name reverseMAthread-package
#' @details
#' This package examines the performance of Mendelian Randomization (MR) and Mediation Analysis methods in the presence of reverse causality.
#' 
#' @section reverseMA for Mendelian Randomization with threading:
#' The reverseMRsim function in the reverseMA R package examines the performance of Mendelian Randomization (MR) methods in the presence of reverse causality. Through simulation studies, this R function examines the type 1 error rate and power for 3 popular MR methods when the role of the intermediate phenotype and outcome were correctly specified and when they were reversed (i.e. reverse causality).
#' 
#' @section Mendelian Randomization Input:
#' nSNP is the number of SNPs generated from a binomial distribution for n subjects (input n) for a given minor allele frequency (input vector MAF).
#' 
#' For the SNPs Xi, the mediator/ exposure M is generated from a normal distribution with the variance (input varM) and the mean as follows:
#' 
#' \eqn{E[M] = \gamma_0 + \Sigma \gamma_x * X_i}
#' 
#' All of these values are inputted by the user (i.e. the intercept gamma0, and the genetic effect size as a vector gammaX).
#' 
#' The outcome Y is generated from a normal distribution with the variance (input varY) and the mean as follows:
#'  
#' \eqn{E[Y] = \beta_0 + \beta_M * M}
#' 
#' All of these values are inputted by the user (i.e. the intercept beta0 and the effect of the mediator directly on the outcome as betaM).
#'
#' After the SNPs X, mediator M, and outcome Y are generated, then the reverseMRsim function compares the power and type 1 error rate of the following 3 methods to detect the path from M to Y: Egger Regression, the Median Weighted Approach, and the Inverse Variance Weighted (IVW) Approach.
#' 
#' @section Mendelian Randomization Example:
#' 
#' For 1,000 subjects (n=1000), we generated 10 SNPs (nSNP=10) with a minor allele frequency of 20\% (specified by MAF) that have a genetic effect size of 0.4 (specified by gammaX) on the normally distributed mediator and the mediator has an effect size varying from 0, 0.2 to 0.3 (specified by betaM) on the normally distributed outcome. We considered 3 MR approaches: Egger Regression, the Median Weighted Approach, and the Inverse Variance Weighted (IVW) Approach.
#' 
#' @section Mendelian Randomization Output:
#' 
#' For the example, we get corresponding plot. In the plot below, the methods ending in NR have the true outcome as the outcome where as the methods ending in R have the true outcome reversed with the mediator. When the mediator and outcome are reversed, the Egger regression and the Median Weighted Approach have an inflated type 1 error rate. While the IVW approach does not have an inflated type 1 error rate, there is very little difference in the IVW approach if the mediator and outcome are reversed, which implies that this approach cannot easily distinguish the causal relationship between the mediator and outcome.
#' 
#' \if{html}{\figure{reverseMRplot.png}{alt="MR_PLOT_image_placeholder"}}
#' \if{latex}{\figure{reverseMRsim.pdf}{alt="MR_PLOT_image_placeholder"}}
#' 
#' 
#' @references  
#' MR.Egger is the Egger Regression approach to MR.
#' \preformatted{
#' Bowden J., Davey Smith G., & Burgess S. (2015). Mendelian Randomization 
#' with invalid instruments: effect estimation and bias detection through 
#' Egger regression. International Journal of Epidemiology, 44(2), 512-525. 
#' }
#' MR.IVW is the Inverse Variant Weighted approach to MR.
#' \preformatted{
#' Burgess, S., Butterworth, A., & Thompson, S. G. (2013). Mendelian 
#' Randomization Analysis With Multiple Genetic Variants Using Summarized 
#' Data. Genetic Epidemiology, 37(7), 658-665.
#' }
#' MR. Median is the Median Weighted approach to MR.
#' \preformatted{
#' Bowden, J., Davey Smith, G., Haycock, P. C., & Burgess, S. (2016). Consistent 
#' Estimation in Mendelian Randomization with Some Invalid Instruments Using a 
#' Weighted Median Estimator. Genetic Epidemiology, 40(4), 304-314. 
#' }
#' 
#' @section reverseMA for Mediation Analysis:
#' These functions examines the performance of mediation analysis methods in the presence of reverse causality.
#' 
#' @section Mediation Analysis Example:
#' 
#' \preformatted{
#' library(reverseMAthread)
#' ?reverseMA # For details on this function
#' 
#' reverseMA(n = 1000, pX = 0.2, gamma0 = 0, gammaX = 0.2, varM = 1, beta0 = 0, betaX = 0, 
#'              betaM = c(0.1, 0.2, 0.3), varY = 1, nSim = 500, nSimImai = 500, SEED = 1, plot.pdf = T, 
#'              plot.name = "reverseMAplot.pdf", alpha_level = 0.05)
#' 
#' reverseMA(n = 1000, pX = 0.2, gamma0 = 0, gammaX = 0, varM = 1, beta0 = 0, betaX = 0.2, 
#'              betaM = c(0.1, 0.2, 0.3), varY = 1, nSim = 500, nSimImai = 500, SEED = 1, plot.pdf = T, 
#'              plot.name = "reverseMAplotDirect.pdf", alpha_level = 0.05)
#' 
#' reverseMA(n = 1000, pX = 0.2, gamma0 = 0, gammaX = 0.2, varM = 1, beta0 = 0, betaX = 0.2, 
#'              betaM = c(0.1, 0.2, 0.3), varY = 1, nSim = 500, nSimImai = 500, SEED = 1, plot.pdf = T, 
#'              plot.name = "reverseMAplotBoth.pdf", alpha_level = 0.05)
#' }
#' @section Mediation Analysis Output:
#' 
#' The example code produces this plot:
#' \if{html}{\figure{reverseMAplot.png}{alt="MA_PLOT_image_placeholder"}}
#' \if{latex}{\figure{reverseMA.pdf}{alt="MA_PLOT_image_placeholder"}}
#' 
#' @section Installation:
#' \preformatted{
#' install.packages("devtools") # devtools must be installed first
#' devtools::install_github("SharonLutz/reverseMAthread")
#' }
NULL
