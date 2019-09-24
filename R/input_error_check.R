input_error_check <- function(n,pX,gamma0,gammaX,varM,beta0,betaX,betaM,varY,
                              nSim,nSimImai,alpha_level, 
                              use_multi_processing, use_cpp, num_jobs){
  
  if(alpha_level>1 | alpha_level<0){stop("Error: alpha_level must be between 0 and 1")}
  if(length(nSimImai) != 1){stop ("Error: nSimImai must be a single integer value")}
  if(nSimImai<=0 || floor(nSimImai)!=ceiling(nSimImai) ){stop("Error: n must be an integer greater than or equal to 1")}
  
  # Error checks
  
  #number of values checks
  if(length(n) != 1){stop ("Error:  n must be a single numeric value")}
  if(length(pX) != 1){stop ("Error: pX must be a single numeric value")}
  if(length(gamma0) != 1){stop ("Error: gamma0 must be a single numeric value")}
  if(length(gammaX) != 1){stop ("Error: gammaX must be a single numeric value")}
  if(length(varM) != 1){stop ("Error: varM must be a single numeric value")}
  if(length(beta0) != 1){stop ("Error: beta0 must be a single numeric value")}
  if(length(betaX) != 1){stop ("Error: betaX must be a single numeric value")}
  if(length(varY) != 1){stop ("Error: varY must be a single numeric value")}
  if(length(nSim) != 1){stop ("Error: nSim must be a single integer value")}
  
  
  if(length(unique(betaM)) != length((betaM))){stop("Error: values in betaM must be unique")}
  if(length(unique(betaM))<2){stop("Error: betaM must be a vector with at least two values")}
  
  #INT style values
  if(n<=0 || floor(n)!=ceiling(n) ){stop("Error: n must be an integer greater than or equal to 1")}
  if(nSim<=0 || floor(nSim)!=ceiling(nSim) ){stop("Error: n must be an integer greater than or equal to 1")}
  
  #valid input checks
  
  if(pX<0 | pX>1){stop("Error: pX must be greater than 0 and less than 1")}
  if(varM<=0){stop("Error: varM must be greater than 0")}
  if(varY<=0){stop("Error: varY must be greater than 0")}
  
  if(use_multi_processing || use_cpp){
    if(use_multi_processing && num_jobs < 2){
      stop("There is no point in using MultiProcessing with less than 2 jobs")
    }
    
    if(parallel::detectCores() == 1){
      warning("your machine may not be suitable for threading or multiprocessing, only 1 core was detected")
    }
    
    if((nSim / num_jobs) < 1.0){
      warning(paste("you don't have enough Simulations in nSim:",nSim," to fully benefit from num_jobs:",num_jobs,sep=""))
    }
  } else {
    if(num_jobs != 1){
      stop("num_jobs != 1 and use_multi_processing=F && use_cpp=F")
    }
  }
}