#' @include simulate_and_mediate.R
#' @import pbapply

mediate_parallel.unix <- function(list_of_job_args, num_jobs=1) {
  options(mc.cores = num_jobs)
  
  result <- pbapply::pblapply(list_of_job_args, simulate_and_mediate, cl = num_jobs)
  attr(result, "dim") <- dim(list_of_job_args)
  
  return(result)
}

mediate_parallel.non_unix <-function(list_of_job_args, num_jobs=1) {
  options(cl.cores = num_jobs)
  snow::setDefaultClusterOptions(type="SOCK")
  this.cluster <- snow::makeCluster(num_jobs)
  on.exit(snow::stopCluster(this.cluster))
  #make sure reverseMAthread is loaded on the nodes
  if(pkgload::is_dev_package("reverseMAthread")){
    #we're in a dev environment, need to load with load_all
    snow::clusterCall(cl=this.cluster,function(){suppressMessages(library(devtools));suppressMessages(load_all())})
  } else {
    #we're being used from an installed copy of reverseMAthread load package explicitly
    snow::clusterCall(cl=this.cluster,function(){suppressMessages(library(reverseMAthread))})
  }
  
  result <- pbapply::pblapply(list_of_job_args, simulate_and_mediate, cl=this.cluster)
  
  dim(result) <- dim(list_of_job_args)
  
  return (result)
}

#' @title Perform Mediation Using Parallel Processing
#' @description starts up multiple processes, configures their environment, and runs simulation and mediation analysis for each item in the input
#' @name mediate_parallel
#' @param list_of_job_args A list or matrix of list objects with named positions that can be used as inputs to simulate_and_mediatie
#' @param num_jobs number of parallel processes to use
#' @return list of \code{MediationProbValues} objects
mediate_parallel <- function(list_of_job_args, num_jobs=1){
  pbapply::pboptions(type="timer", style=1)
  if(.Platform$OS.type == "unix") {
    result <- mediate_parallel.unix(list_of_job_args=list_of_job_args, num_jobs=num_jobs)
  } else {
    result <- mediate_parallel.non_unix(list_of_job_args=list_of_job_args, num_jobs=num_jobs)
  }
  return(result)
}
