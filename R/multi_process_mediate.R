source("R/mediate_variables_class.R")
suppressMessages(library(parallel))

simulate_and_mediate <- function(med_model_vars){
  suppressMessages(library(mediation))
  #my test begins here
  
  models <- assembleLinearModels(med_model_vars)
  
  # Fit the mediation model
  set.seed(med_model_vars@SEED)
  
  med.out <- mediate(models@med.fit, models@out.fit, treat = "X",mediator = "M1",sims = nSimImai)
  med.out.r <- mediate(models@med.fit.r, models@out.fit.r, treat = "X",mediator = "M2",sims = nSimImai)
  
  summary_obj = summary(med.out)
  summary_obj.r = summary(med.out.r)
  
  return(MediationProbValues(
    pval_direct = summary_obj$z.avg.p, 
    pval_indirect = summary_obj$d.avg.p, 
    pval_direct_r = summary_obj.r$z.avg.p,
    pval_indirect_r = summary_obj.r$d.avg.p)
    )
}

if (.Platform$OS.type == "unix") {
  library(parallel)
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
    result <- mclapply(list_of_job_args, simulate_and_mediate, mc.silent = TRUE)
    attr(result, "dim") <- dim(list_of_job_args)
    if(is.null(old_nSimImai)){
      rm("nSimImai", envir=g_env)
    } else {
      g_env[["nSimImai"]] = old_nSimImai
    }
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