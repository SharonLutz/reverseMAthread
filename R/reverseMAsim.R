#' @include reverseMAsim_impl.R

#' @export
#' @title reverseMAsim
#' @description A function to simulate the performance of the mediate function from the mediation package in scenarios of reverse causality.
#' @author Annie Thwing, Sharon Lutz
#' @param n is the sample size.
#' @param pX is the minor allele frequency
#' @param gamma0 is the intercept for M
#' @param gammaX is the association of X with M
#' @param varM is the variance of M
#' @param beta0 is the intercept for Y
#' @param betaX is the direct effect of X on Y
#' @param betaM is a vector of different associations of M with Y
#' @param varY is the variance of Y
#' @param nSim is the number of simulations to run
#' @param nSimImai is the number of simulations to run in mediate from the mediation package
#' @param SEED is the seed
#' @param plot.pdf is T to output a plot, is F to not output a plot
#' @param plot.name is the name of the plot
#' @param alpha_level is the significance level
#' @param use_multi_processing use multi_processing instead of threading for a speed boost
#' @param use_cpp use Rcpp with Eigen
#' @param num_jobs the number of cores to use, i.e. the number of parallel procesess or threads to spawn
#' @return a matrix of the power of the mediate method from the mediation package to detect an effect of the mediator M on the outcome Y when M and Y are correctly specified and also when they are incorrectly specified (the true mediator is Y and the true outcome is M)
reverseMAsim <-
function(n=1000,pX=0.2,gamma0=0,gammaX=0.1,varM=1,beta0=0,betaX=1,betaM=c(0,0.1,0.2),varY=1,
         nSim=100,nSimImai=1000,SEED=1,plot.pdf=T,plot.name="reverseMAsim.pdf",alpha_level=0.05, 
         use_multi_processing=F, use_cpp=F, num_jobs=1){
  # cat(paste("num_jobs:",num_jobs, "; use_cpp:",use_cpp, "; use_multi_processing:",use_multi_processing, "\n",sep = ""))
  pbapply::pboptions(style=1,type="timer")
  if(use_multi_processing){
    reverseMAsim.MultiProcess(n=n, pX=pX,gamma0=gamma0, gammaX=gammaX, varM=varM, beta0=beta0, betaX=betaX,
                              betaM=betaM, varY=varY, nSim=nSim, nSimImai=nSimImai, SEED=SEED, plot.pdf=plot.pdf,
                              plot.name=plot.name, alpha_level=alpha_level,use_cpp=use_cpp,
                              num_jobs=num_jobs)
  } else {
    reverseMAsim.SingleProcess(n=n, pX=pX,gamma0=gamma0, gammaX=gammaX, varM=varM, beta0=beta0, betaX=betaX,
                               betaM=betaM, varY=varY, nSim=nSim, nSimImai=nSimImai, SEED=SEED, plot.pdf=plot.pdf,
                               plot.name=plot.name, alpha_level=alpha_level,use_cpp=use_cpp,
                               num_jobs=num_jobs)
  }
}
