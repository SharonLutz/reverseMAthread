#' @include mediate_s4_classes.R
#' @include multi_process_mediate_impl.R
#' @importFrom mediation mediate

#' @title Perform Mediation Using Parallel Processing
#' @description starts up multiple processes, configures their environment, and runs simulation and mediation analysis for each item in the input
#' @name mediate_parallel
#' @param list_of_job_args A list of \code{MediateModelVariables} objects
#' @param nSimImai number of mediation sims to perform
#' @return list of \code{MediationProbValues} objects
mediate_parallel <- function(list_of_job_args, nSimImai=1000, use_cpp=F, num_jobs=getOption("mediate.jobs", parallel::detectCores() - 1)){
  pbapply::pboptions(type="timer", style=1)
  if(.Platform$OS.type == "unix") {
    result <- mediate_parallel.unix(list_of_job_args=list_of_job_args, nSimImai=nSimImai, use_cpp=use_cpp, num_jobs=num_jobs)
  } else {
    result <- mediate_parallel.non_unix(list_of_job_args=list_of_job_args, nSimImai=nSimImai, use_cpp=use_cpp, num_jobs=num_jobs)
  }
  return(result)
}
