source("mediate_variables_class.R")

simulate_and_mediate <- function(med_vars = NULL){
  library(mediation)
  #my test begins here
  
  if(is.null(med_vars)){
    med_vars <- MediateVariables()
  }
  
  data_vars <- generateData(med_vars)
  
  models <- assembleLinearModels(data_vars)
  
  # Fit the mediation model
  
  med.out <- mediate(models@med.fit, models@out.fit, treat = "X",mediator = "M",sims = med_vars@nSimImai)
  
  return(med.out)
}

if (.Platform$OS.type == "unix") {
  library(parallel)
  do_mediate <- function(list_of_job_args){
    options(mc.cores = getOption("mediate.cores", detectCores() - 1))
    result <- mclapply(list_of_job_args, simulate_and_mediate, mc.silent = TRUE)
    return(result)
  }
}else{
  library(parallel)
  library(snow)
  
  #wrap input args somehow into a list to supply to cluster nodes
  #ensure the function called will return all required information wrapped up so it will fit in a list
  #maybe could use environment objects....
  do_mediate <- function(list_of_job_args){
    num_cores <- getOption("mediate.cores", detectCores() - 1)
    options(cl.cores = num_cores)
    this.cluster <- makeCluster(num_cores)
    clusterCall(cl=this.cluster,function(){source("mediate_variables_class.R")})
    
    result <- parLapply(cl=this.cluster, list_of_job_args, simulate_and_mediate)
    on.exit(stopCluster(this.cluster))
    
    return (result)
  }
  
}