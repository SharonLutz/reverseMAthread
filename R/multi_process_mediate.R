source("R/mediate_variables_class.R")
suppressMessages(library(parallel))

simulate_and_mediate <- function(med_data = NULL){
  suppressMessages(library(mediation))
  #my test begins here
  
  if(is.null(med_data)){
    med_vars = MediateVariables()
    med_data = generateData(med_vars,1)
  }
  
  models <- assembleLinearModels(med_data)
  
  # Fit the mediation model
  
  med.out <- mediate(models@med.fit, models@out.fit, treat = "X",mediator = "M",sims = nSimImai)
  summary_obj = summary(med.out)
  return(MediationProbValues(pval_direct=summary_obj$z.avg.p, pval_indirect=summary_obj$d.avg.p))
}

if (.Platform$OS.type == "unix") {
  library(parallel)
  mediate_parallel <- function(list_of_job_args, nSimImai=1000){
    options(mc.cores = getOption("mediate.cores", detectCores() - 1))
    # TODO: make sure nSimImai value exists for all processes created by mclapply
    result <- mclapply(list_of_job_args, simulate_and_mediate, mc.silent = TRUE)
    attr(result, "dim") <- dim(list_of_job_args)
    return(result)
  }
}else{
  suppressMessages(library(snow))
  
  #wrap input args somehow into a list to supply to cluster nodes
  #ensure the function called will return all required information wrapped up so it will fit in a list
  #maybe could use environment objects....
  mediate_parallel <- function(list_of_job_args, nSimImai=1000){
    num_cores <- getOption("mediate.cores", detectCores() - 1)
    options(cl.cores = num_cores)
    this.cluster <- makeCluster(num_cores)
    on.exit(stopCluster(this.cluster))
    clusterCall(cl=this.cluster,function(){suppressMessages(source("R/mediate_variables_class.R"))})
    clusterExport(cl=this.cluster,c("nSimImai"),envir=environment())
    
    result <- parLapply(cl=this.cluster, list_of_job_args, simulate_and_mediate)
    attr(result, "dim") <- dim(list_of_job_args)
    return (result)
  }
  
}