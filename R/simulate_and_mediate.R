
#' @title Simulate Data and Mediate Linear Models
#' @description runs forward and reverse mediation on argument contents
#' @param data_element list with named elements that will be used to assemble linear models to use with mediation
#' @return a \code{MediationProbValues} instance with values obtained by running mediate on forward and reverse linear models.
simulate_and_mediate <- function(data_element){
  use_cpp = data_element[["use_cpp"]]
  num_jobs = data_element[["num_jobs"]]
  use_multi_processing = data_element[["use_multi_processing"]]
  
  nSimImai = data_element[["nSimImai"]]
  
  df = data.frame(X=data_element[["X"]], Y1=data_element[["Y"]], M1=data_element[["M"]], M2=data_element[["Y"]], Y2=data_element[["M"]])
  
  med.fit = lm("M1~X", data=df)
  out.fit = lm("Y1~X+M1", data=df)
  med.fit.r = lm("M2~X", data=df)
  out.fit.r = lm("Y2~X+M2", data=df)
  
  if(exists(".Random.seed",envir = .GlobalEnv) && !is.null(.GlobalEnv[[".Random.seed"]])){
    old_rand_state = .GlobalEnv[[".Random.seed"]]
  }
  
  .GlobalEnv[[".Random.seed"]] = data_element[["RAND_STATE"]]
  
  if(use_cpp){
    
    if(use_multi_processing){
      options(mediate.threads = 1)
    } else {
      options(mediate.threads = num_jobs)
    }
    
    med.out <- mediate_with_rcpp(med.fit, out.fit, treat = "X",mediator = "M1",sims = nSimImai)
    med.out.r <- mediate_with_rcpp(med.fit.r, out.fit.r, treat = "X",mediator = "M2",sims = nSimImai)
  } else {
    med.out <- mediation::mediate(med.fit, out.fit, treat = "X",mediator = "M1",sims = nSimImai)
    med.out.r <- mediation::mediate(med.fit.r, out.fit.r, treat = "X",mediator = "M2",sims = nSimImai)
  }
  
  summary_obj = mediation::summary.mediate(med.out)
  summary_obj.r = mediation::summary.mediate(med.out.r)
  
  result = list(pval_direct = summary_obj[["z.avg.p"]],
                pval_indirect = summary_obj[["d.avg.p"]],
                pval_direct_r = summary_obj.r[["z.avg.p"]],
                pval_indirect_r = summary_obj.r[["d.avg.p"]])
  
  if(!is.null(old_rand_state)){
    .GlobalEnv[[".Random.seed"]] = old_rand_state
  }
  
  return(result)
}