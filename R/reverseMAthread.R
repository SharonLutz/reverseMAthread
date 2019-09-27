#' @include reverseMAthread_impl.R

#' @export
reverseMAthread <-
  function(n=1000,pX=0.2,gamma0=0,gammaX=0.1,varM=1,beta0=0,betaX=1,betaM=c(0,0.1,0.2),varY=1,
           nSim=100,nSimImai=1000,SEED=1,plot.pdf=T,plot.name="reverseMAthread.pdf",alpha_level=0.05, 
           use_multi_processing=F, use_cpp=F, num_jobs=1){
    
    # cat(paste("num_jobs:",num_jobs, "; use_cpp:",use_cpp, "; use_multi_processing:",use_multi_processing, "\n",sep = ""))
    pbapply::pboptions(style=1,type="timer")
    if(use_multi_processing){
      reverseMAthread.MultiProcess(n=n, pX=pX,gamma0=gamma0, gammaX=gammaX, varM=varM, beta0=beta0, betaX=betaX,
                                betaM=betaM, varY=varY, nSim=nSim, nSimImai=nSimImai, SEED=SEED, plot.pdf=plot.pdf,
                                plot.name=plot.name, alpha_level=alpha_level,use_cpp=use_cpp,
                                num_jobs=num_jobs)
    } else {
      reverseMAthread.SingleProcess(n=n, pX=pX,gamma0=gamma0, gammaX=gammaX, varM=varM, beta0=beta0, betaX=betaX,
                                 betaM=betaM, varY=varY, nSim=nSim, nSimImai=nSimImai, SEED=SEED, plot.pdf=plot.pdf,
                                 plot.name=plot.name, alpha_level=alpha_level,use_cpp=use_cpp,
                                 num_jobs=num_jobs)
    }
  }