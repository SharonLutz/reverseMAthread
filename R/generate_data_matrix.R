generate_data_matrix <- function(
  n=1000,pX=0.2,
  gamma0=0,gammaX=0.1,
  varM=1,varY=1,
  beta0=0,betaX=1,betaM=c(0,0.1,0.2),
  nSim=100, nSimImai=1000, SEED=1,
  use_multi_processing=F, use_cpp=F, num_jobs=1
  ){
  
  # Error checks
  if(n<0 | n==0 | floor(n)!=ceiling(n) ){stop("Error: n must be an integer greater than or equal to 1")}
  if(pX<0 | pX>1){stop("Error: pX must be greater than 0 and less than 1")}
  if(!varM>0){stop("Error: varM must be greater than 0")}
  if(!varY>0){stop("Error: varY must be greater than 0")}
  if(length(unique(betaM))<2){stop("Error: betaM must be a vector with at least two values")}
  if(length(nSimImai) != 1){stop ("Error: nSimImai must be a single integer value")}
  if(nSimImai<=0 || floor(nSimImai)!=ceiling(nSimImai) ){stop("Error: n must be an integer greater than or equal to 1")}
  
  num_betaM = length(betaM)
  
  mat_total <- matrix(0,nrow=num_betaM,ncol=4)
  colnames(mat_total) <- c("DirectNR","IndirectNR","DirectR","IndirectR")
  
  #generate the data needed to make linear models.
  data_matrix = matrix(list(), nrow=num_betaM, ncol=nSim)
  
  sd_M = sqrt(varM)
  sd_Y = sqrt(varY)
  old_rand_state = NULL
  
  if(exists(".Random.seed",envir = .GlobalEnv) && !is.null(.GlobalEnv[[".Random.seed"]])){
    old_rand_state = .GlobalEnv[[".Random.seed"]]
  }
  
  next_seed = SEED
  
  for(i in 1:nSim){
    
    for(bM.ind in 1:num_betaM){
      set.seed(next_seed)
      
      X = rbinom(n, 2, pX)
      X[X==2]<-1 #dominant model
      M = rnorm(n=n, mean = (gamma0 + gammaX * X), sd=sd_M)
      Y = rnorm(n=n, mean = (beta0 + betaX * X + betaM[[bM.ind]] * M), sd=sd_Y)
      
      data_element = list( X = X,
                           M = M,
                           Y = Y,
                           nSimImai=nSimImai,
                           RAND_STATE = .GlobalEnv[[".Random.seed"]],
                           use_cpp=use_cpp,
                           use_multi_processing=use_multi_processing,
                           num_jobs=num_jobs
                           )
      
      next_seed = next_seed + 1
      
      data_matrix[[bM.ind, i]] = data_element
    }
    
  }
  
  if(!is.null(old_rand_state)){
    .GlobalEnv[[".Random.seed"]] = old_rand_state
  }
  
  return(data_matrix)
  
}