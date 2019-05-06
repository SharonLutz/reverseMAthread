#' @include mediate_s4_classes.R
#' @importFrom mediation mediate



#' @title Simulate Data and Mediate Linear Models
#' @description runs forward and reverse mediation on argument contents
#' @param med_model_vars an instance of \code{MediateModelVariables}
#' @return a \code{MediationProbValues} instance with values obtained by running mediate on forward and reverse linear models.
simulate_and_mediate <- function(med_model_vars){
  #my test begins here
  
  models <- assembleLinearModels(med_model_vars)
  
  # Fit the mediation model
  set.seed(med_model_vars@SEED)
  
  med.out <- stripped_down_mediate_with_rcpp(models@med.fit, models@out.fit, treat = "X",mediator = "M1",sims = nSimImai)
  med.out.r <- stripped_down_mediate_with_rcpp(models@med.fit.r, models@out.fit.r, treat = "X",mediator = "M2",sims = nSimImai)
  
  # summary_obj = summary(med.out)
  # summary_obj.r = summary(med.out.r)
  # 
  # return(MediationProbValues(
  #   pval_direct = summary_obj$z.avg.p, 
  #   pval_indirect = summary_obj$d.avg.p, 
  #   pval_direct_r = summary_obj.r$z.avg.p,
  #   pval_indirect_r = summary_obj.r$d.avg.p)
  #   )
  return(
    MediationProbValues(
      pval_direct = med.out@direct_p,
      pval_indirect = med.out@indirect_p,
      pval_direct_r = med.out.r@direct_p,
      pval_indirect_r = med.out.r@indirect_p
    )
  )
}

#' @title Perform Mediation Using Parallel Processing
#' @description starts up multiple processes, configures their environment, and runs simulation and mediation analysis for each item in the input
#' @name mediate_parallel
#' @param list_of_job_args A list of \code{MediateModelVariables} objects
#' @param nSimImai number of mediation sims to perform
#' @return list of \code{MediationProbValues} objects
if (.Platform$OS.type == "unix") {
  mediate_parallel <- function(list_of_job_args, nSimImai=1000){
    g_env = globalenv()
    if(exists("nSimImai", envir = g_env)){
      old_nSimImai = g_env[["nSimImai"]]
    } else {
      old_nSimImai = NULL
    }
    g_env[["nSimImai"]]=nSimImai
    options(mc.cores = getOption("mediate.cores", detectCores() - 1))
    # TODO: make sure nSimImai value exists for all processes created by mclapply
    result <- parallel::mclapply(list_of_job_args, simulate_and_mediate, mc.silent = TRUE)
    attr(result, "dim") <- dim(list_of_job_args)
    if(is.null(old_nSimImai)){
      rm("nSimImai", envir=g_env)
    } else {
      g_env[["nSimImai"]] = old_nSimImai
    }
    return(result)
  }
}else{
  #wrap input args somehow into a list to supply to cluster nodes
  #ensure the function called will return all required information wrapped up so it will fit in a list
  #maybe could use environment objects....
  mediate_parallel <- function(list_of_job_args, nSimImai=1000){
    num_cores <- getOption("mediate.cores", detectCores() - 1)
    options(cl.cores = num_cores)
    this.cluster <- snow::makeCluster(num_cores)
    on.exit(snow::stopCluster(this.cluster))
    #make sure reverseC is loaded on the nodes
    if(pkgload::is_dev_package("reverseC")){
      #we're in a dev environment, need to load with load_all
      snow::clusterCall(cl=this.cluster,function(){suppressMessages(library(devtools));suppressMessages(load_all())})
    } else {
      #we're being used from an installed copy of reverseC load package explicitly
      snow::clusterCall(cl=this.cluster,function(){suppressMessages(library(reverseC))})
    }
    
    snow::clusterExport(cl=this.cluster,c("nSimImai"),envir=environment())
    
    result <- snow::parLapply(cl=this.cluster, list_of_job_args, simulate_and_mediate)
    
    attr(result, "dim") <- dim(list_of_job_args)
    
    return (result)
  }
}