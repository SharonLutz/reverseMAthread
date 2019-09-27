#' @include mediate_s4_classes.R
#' @include multi_process_mediate_impl.R

mediate_parallel <- function(list_of_job_args, nSimImai=1000, use_cpp=F, num_jobs=getOption("mediate.jobs", parallel::detectCores() - 1)){
  pbapply::pboptions(type="timer", style=1)
  if(.Platform$OS.type == "unix") {
    result <- mediate_parallel.unix(list_of_job_args=list_of_job_args, nSimImai=nSimImai, use_cpp=use_cpp, num_jobs=num_jobs)
  } else {
    result <- mediate_parallel.non_unix(list_of_job_args=list_of_job_args, nSimImai=nSimImai, use_cpp=use_cpp, num_jobs=num_jobs)
  }
  return(result)
}
