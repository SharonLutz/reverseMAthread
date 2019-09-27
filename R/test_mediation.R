#our calls always this style: mediate(model.m, model.y, treat= "X", mediator="M", sims = nSimImai)
#both models are lm with no special types / conditions
test.mediate.vanilla <- 
  function(model.m, model.y, sims = 1000, treat = "treat.name", mediator = "med.name",
           conf.level = .95, control.value = 0, treat.value = 1, return_context=F, context_before=F, export_loop_vars=F){
    
    # Model frames for M and Y models
    m.data <- model.frame(model.m)  # Call.M$data
    y.data <- model.frame(model.y)  # Call.Y$data
    
    # Specify group names
    group.m <- NULL
    group.y <- NULL
    group.out <- NULL
    group.id.m <- NULL
    group.id.y <- NULL
    group.id <- NULL
    group.name <- NULL
    
    # Numbers of observations and categories
    n.m <- nrow(m.data)
    n.y <- nrow(y.data)
    
    if(n.m != n.y){
      stop("number of observations do not match between mediator and outcome models")
    } else{
      n <- n.m
    }
    m <- length(sort(unique(model.frame(model.m)[,1])))
    
    
    # Extracting weights from models
    weights.m <- model.weights(m.data)
    weights.y <- model.weights(y.data)
    weights.m <- rep(1,nrow(m.data))
    weights.y <- rep(1,nrow(y.data))
    weights <- weights.m
    
    cat.0 <- control.value
    cat.1 <- treat.value
    
    ########################################################################
    ## Case I-1: Quasi-Bayesian Monte Carlo
    ########################################################################
    
    # Get mean and variance parameters for mediator simulations
    MModel.coef <- coef(model.m)
    scalesim.m <- FALSE
    
    MModel.var.cov <- vcov(model.m)
    
    YModel.coef <- coef(model.y)
    scalesim.y <- FALSE
    
    YModel.var.cov <- vcov(model.y)
    
    if(sum(is.na(MModel.coef)) > 0){
      stop("NA in model coefficients; rerun models with nonsingular design matrix")
    }
    MModel <- mvtnorm::rmvnorm(sims, mean=MModel.coef, sigma=MModel.var.cov)
    
    if(sum(is.na(YModel.coef)) > 0){
      stop("NA in model coefficients; rerun models with nonsingular design matrix")
    }
    YModel <- mvtnorm::rmvnorm(sims, mean=YModel.coef, sigma=YModel.var.cov)
    
    #####################################
    ##  Mediator Predictions
    #####################################
    
    pred.data.t <- pred.data.c <- m.data
    
    pred.data.t[,treat] <- cat.1
    pred.data.c[,treat] <- cat.0
    
    mmat.t <- model.matrix(terms(model.m), data=pred.data.t)
    mmat.c <- model.matrix(terms(model.m), data=pred.data.c)
    
    sigma <- summary(model.m)$sigma
    error <- rnorm(sims*n, mean=0, sd=sigma)
    muM1 <- tcrossprod(MModel, mmat.t)
    muM0 <- tcrossprod(MModel, mmat.c)
    PredictM1 <- muM1 + matrix(error, nrow=sims)
    PredictM0 <- muM0 + matrix(error, nrow=sims)
    
    rm(error)
    rm(mmat.t, mmat.c)
    
    #####################################
    ##  Outcome Predictions
    #####################################
    if(return_context && context_before){
      return(environment())
    }
    effects.tmp <- array(NA, dim = c(n, sims, 4))
    
    exported = F
    
    for(e in 1:4){
      tt <- switch(e, c(1,1,1,0), c(0,0,1,0), c(1,0,1,1), c(1,0,0,0))
      Pr1 <- matrix(nrow=n, ncol=sims)
      Pr0 <- matrix(nrow=n, ncol=sims)
      
      for(j in 1:sims){
        pred.data.t <- pred.data.c <- y.data
        
        # Set treatment values
        cat.t <- ifelse(tt[1], cat.1, cat.0)
        cat.c <- ifelse(tt[2], cat.1, cat.0)
        cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
        cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)
        
        pred.data.t[,treat] <- cat.t
        pred.data.c[,treat] <- cat.c
        
        # Set mediator values
        PredictMt <- PredictM1[j,] * tt[3] + PredictM0[j,] * (1 - tt[3])
        PredictMc <- PredictM1[j,] * tt[4] + PredictM0[j,] * (1 - tt[4])
        
        pred.data.t[,mediator] <- PredictMt
        pred.data.c[,mediator] <- PredictMc
        
        ymat.t <- model.matrix(terms(model.y), data=pred.data.t)
        ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
        
        # ymodel_j = YModel[j,]
        # ymodel_j_names = list()
        # ymodel_j_names[[1]] = NULL
        # ymodel_j_names[[2]] = names(ymodel_j)
        # ymodel_j_mat = matrix(ymodel_j, nrow=1, dimnames=ymodel_j_names)
        # 
        # Pr1[,j] <- ymodel_j_mat %*% t(ymat.t)
        # Pr0[,j] <- ymodel_j_mat %*% t(ymat.c)
        
        Pr1[,j] <- t(as.matrix(YModel[j,])) %*% t(ymat.t)
        Pr0[,j] <- t(as.matrix(YModel[j,])) %*% t(ymat.c)
        if((!exported) && export_loop_vars){
          genv = globalenv()
          genv[["vanillaR.ymat.t"]] = ymat.t
          genv[["vanillaR.ymat.c"]] = ymat.c
          genv[["vanillaR.pred.data.t"]] = pred.data.t
          genv[["vanillaR.pred.data.c"]] = pred.data.c
          genv[["vanillaR.tt"]] = tt
          genv[["vanillaR.PredictMt"]] = PredictMt
          genv[["vanillaR.PredictMc"]] = PredictMc
          genv[["vanillaR.cat.t"]] = cat.t
          genv[["vanillaR.cat.c"]] = cat.c
          genv[["vanillaR.cat.t.ctrl"]] = cat.t.ctrl
          genv[["vanillaR.cat.c.ctrl"]] = cat.c.ctrl
          exported= T
        }
        
        rm(ymat.t, ymat.c, pred.data.t, pred.data.c)
      }
      
      effects.tmp[,,e] <- Pr1 - Pr0 ### e=1:mediation(1); e=2:mediation(0); e=3:direct(1); e=4:direct(0)
      if(export_loop_vars){
        genv = globalenv()
        genv[[paste("vanillaR.Pr1", e ,sep=".")]] = Pr1
        genv[[paste("vanillaR.Pr0", e ,sep=".")]] = Pr0
      }
      rm(Pr1, Pr0)
    }
    
    rm(PredictM1, PredictM0, YModel, MModel)
    
    et1<-effects.tmp[,,1] ### mediation effect (1)
    et2<-effects.tmp[,,2] ### mediation effect (0)
    et3<-effects.tmp[,,3] ### direct effect (1)
    et4<-effects.tmp[,,4] ### direct effect (0)
    
    if(return_context && !context_before){
      return(environment())
    }
    
    delta.1 <- t(as.matrix(apply(et1, 2, weighted.mean, w=weights)))
    delta.0 <- t(as.matrix(apply(et2, 2, weighted.mean, w=weights)))
    zeta.1 <- t(as.matrix(apply(et3, 2, weighted.mean, w=weights)))
    zeta.0 <- t(as.matrix(apply(et4, 2, weighted.mean, w=weights)))
    rm(effects.tmp)
    
    tau <- (zeta.1 + delta.0 + zeta.0 + delta.1)/2
    nu.0 <- delta.0/tau
    nu.1 <- delta.1/tau
    delta.avg <- (delta.1 + delta.0)/2
    zeta.avg <- (zeta.1 + zeta.0)/2
    nu.avg <- (nu.1 + nu.0)/2
    
    d0 <- mean(delta.0)			# mediation effect
    d1 <- mean(delta.1)
    z1 <- mean(zeta.1)			# direct effect
    z0 <- mean(zeta.0)
    tau.coef <- mean(tau)	  	        # total effect
    n0 <- median(nu.0)
    n1 <- median(nu.1)
    d.avg <- (d0 + d1)/2
    z.avg <- (z0 + z1)/2
    n.avg <- (n0 + n1)/2
    
    ########################################################################
    ## Compute Outputs and Put Them Together
    ########################################################################
    
    low <- (1 - conf.level)/2
    high <- 1 - low
    
    d0.ci <- quantile(delta.0, c(low,high), na.rm=TRUE)
    d1.ci <- quantile(delta.1, c(low,high), na.rm=TRUE)
    tau.ci <- quantile(tau, c(low,high), na.rm=TRUE)
    z1.ci <- quantile(zeta.1, c(low,high), na.rm=TRUE)
    z0.ci <- quantile(zeta.0, c(low,high), na.rm=TRUE)
    n0.ci <- quantile(nu.0, c(low,high), na.rm=TRUE)
    n1.ci <- quantile(nu.1, c(low,high), na.rm=TRUE)
    d.avg.ci <- quantile(delta.avg, c(low,high), na.rm=TRUE)
    z.avg.ci <- quantile(zeta.avg, c(low,high), na.rm=TRUE)
    n.avg.ci <- quantile(nu.avg, c(low,high), na.rm=TRUE)
    
    # p-values
    d0.p <- mediation:::pval(delta.0, d0)
    d1.p <- mediation:::pval(delta.1, d1)
    d.avg.p <- mediation:::pval(delta.avg, d.avg)
    z0.p <- mediation:::pval(zeta.0, z0)
    z1.p <- mediation:::pval(zeta.1, z1)
    z.avg.p <- mediation:::pval(zeta.avg, z.avg)        
    n0.p <- mediation:::pval(nu.0, n0)
    n1.p <- mediation:::pval(nu.1, n1)
    n.avg.p <- mediation:::pval(nu.avg, n.avg)
    tau.p <- mediation:::pval(tau, tau.coef)
    
    # Detect whether models include T-M interaction
    INT <- paste(treat,mediator,sep=":") %in% attr(terms(model.y),"term.labels") |
      paste(mediator,treat,sep=":") %in% attr(terms(model.y),"term.labels")
    
    return(SimpleMediateResult(direct_p = z.avg.p, indirect_p=d.avg.p))
  }

#our calls always this style: mediate(model.m, model.y, treat= "X", mediator="M", sims = nSimImai)
#both models are lm with no special types / conditions
test.mediate.rcpp <- 
  function(model.m, model.y, sims = 1000, treat = "treat.name", mediator = "med.name",
           conf.level = .95, control.value = 0, treat.value = 1, num_cores=1, return_context=F, context_before=F, export_loop_vars=F){
    if(export_loop_vars && (num_cores > 1)){
      stop("Cannot export loop vars if using more than num_cores==1")
    }
    # Model frames for M and Y models
    m.data <- model.frame(model.m)  # Call.M$data
    y.data <- model.frame(model.y)  # Call.Y$data
    
    # Specify group names
    group.m <- NULL
    group.y <- NULL
    group.out <- NULL
    group.id.m <- NULL
    group.id.y <- NULL
    group.id <- NULL
    group.name <- NULL
    
    # Numbers of observations and categories
    n.m <- nrow(m.data)
    n.y <- nrow(y.data)
    
    if(n.m != n.y){
      stop("number of observations do not match between mediator and outcome models")
    } else{
      n <- n.m
    }
    m <- length(sort(unique(model.frame(model.m)[,1])))
    
    
    # Extracting weights from models
    weights.m <- model.weights(m.data)
    weights.y <- model.weights(y.data)
    weights.m <- rep(1,nrow(m.data))
    weights.y <- rep(1,nrow(y.data))
    weights <- weights.m
    
    cat.0 <- control.value
    cat.1 <- treat.value
    
    ########################################################################
    ## Case I-1: Quasi-Bayesian Monte Carlo
    ########################################################################
    
    # Get mean and variance parameters for mediator simulations
    MModel.coef <- coef(model.m)
    scalesim.m <- FALSE
    
    MModel.var.cov <- vcov(model.m)
    
    YModel.coef <- coef(model.y)
    scalesim.y <- FALSE
    
    YModel.var.cov <- vcov(model.y)
    
    if(sum(is.na(MModel.coef)) > 0){
      stop("NA in model coefficients; rerun models with nonsingular design matrix")
    }
    MModel <- mvtnorm::rmvnorm(sims, mean=MModel.coef, sigma=MModel.var.cov)
    
    if(sum(is.na(YModel.coef)) > 0){
      stop("NA in model coefficients; rerun models with nonsingular design matrix")
    }
    YModel <- mvtnorm::rmvnorm(sims, mean=YModel.coef, sigma=YModel.var.cov)
    
    #####################################
    ##  Mediator Predictions
    #####################################
    
    pred.data.t <- pred.data.c <- m.data
    
    pred.data.t[,treat] <- cat.1
    pred.data.c[,treat] <- cat.0
    
    mmat.t <- model.matrix(terms(model.m), data=pred.data.t)
    mmat.c <- model.matrix(terms(model.m), data=pred.data.c)
    
    sigma <- summary(model.m)$sigma
    error <- rnorm(sims*n, mean=0, sd=sigma)
    muM1 <- tcrossprod(MModel, mmat.t)
    muM0 <- tcrossprod(MModel, mmat.c)
    PredictM1 <- muM1 + matrix(error, nrow=sims)
    PredictM0 <- muM0 + matrix(error, nrow=sims)
    
    rm(error)
    rm(mmat.t, mmat.c)
    
    #####################################
    ##  Outcome Predictions
    #####################################
    
    if(return_context && context_before){
      return(environment())
    }
    
    if(num_cores > 1){
      threaded_mediate_helper(environment(), num_cores)
    } else {
      if(export_loop_vars){
        mediate_helper_variable_exporter(environment())
      } else {
        mediate_helper(environment())
      }
    }
    
    if(return_context && !context_before){
      return(environment())
    }
    
    delta.1 <- t(as.matrix(apply(et1, 2, weighted.mean, w=weights)))
    delta.0 <- t(as.matrix(apply(et2, 2, weighted.mean, w=weights)))
    zeta.1 <- t(as.matrix(apply(et3, 2, weighted.mean, w=weights)))
    zeta.0 <- t(as.matrix(apply(et4, 2, weighted.mean, w=weights)))
    
    tau <- (zeta.1 + delta.0 + zeta.0 + delta.1)/2
    nu.0 <- delta.0/tau
    nu.1 <- delta.1/tau
    delta.avg <- (delta.1 + delta.0)/2
    zeta.avg <- (zeta.1 + zeta.0)/2
    nu.avg <- (nu.1 + nu.0)/2
    
    d0 <- mean(delta.0)			# mediation effect
    d1 <- mean(delta.1)
    z1 <- mean(zeta.1)			# direct effect
    z0 <- mean(zeta.0)
    tau.coef <- mean(tau)	  	        # total effect
    n0 <- median(nu.0)
    n1 <- median(nu.1)
    d.avg <- (d0 + d1)/2
    z.avg <- (z0 + z1)/2
    n.avg <- (n0 + n1)/2
    
    ########################################################################
    ## Compute Outputs and Put Them Together
    ########################################################################
    
    low <- (1 - conf.level)/2
    high <- 1 - low
    
    d0.ci <- quantile(delta.0, c(low,high), na.rm=TRUE)
    d1.ci <- quantile(delta.1, c(low,high), na.rm=TRUE)
    tau.ci <- quantile(tau, c(low,high), na.rm=TRUE)
    z1.ci <- quantile(zeta.1, c(low,high), na.rm=TRUE)
    z0.ci <- quantile(zeta.0, c(low,high), na.rm=TRUE)
    n0.ci <- quantile(nu.0, c(low,high), na.rm=TRUE)
    n1.ci <- quantile(nu.1, c(low,high), na.rm=TRUE)
    d.avg.ci <- quantile(delta.avg, c(low,high), na.rm=TRUE)
    z.avg.ci <- quantile(zeta.avg, c(low,high), na.rm=TRUE)
    n.avg.ci <- quantile(nu.avg, c(low,high), na.rm=TRUE)
    
    # p-values
    d0.p <- mediation:::pval(delta.0, d0)
    d1.p <- mediation:::pval(delta.1, d1)
    d.avg.p <- mediation:::pval(delta.avg, d.avg)
    z0.p <- mediation:::pval(zeta.0, z0)
    z1.p <- mediation:::pval(zeta.1, z1)
    z.avg.p <- mediation:::pval(zeta.avg, z.avg)        
    n0.p <- mediation:::pval(nu.0, n0)
    n1.p <- mediation:::pval(nu.1, n1)
    n.avg.p <- mediation:::pval(nu.avg, n.avg)
    tau.p <- mediation:::pval(tau, tau.coef)
    
    # Detect whether models include T-M interaction
    INT <- paste(treat,mediator,sep=":") %in% attr(terms(model.y),"term.labels") |
      paste(mediator,treat,sep=":") %in% attr(terms(model.y),"term.labels")
    
    return(SimpleMediateResult(direct_p = z.avg.p, indirect_p=d.avg.p))
  }

export_environment <- function(env){
  glob_env = globalenv()
  for(item in names(env)){
    glob_env[[item]] = env[[item]]
  }
}

collect_vanillaR_loop_exports <- function(){
  local_env = new.env(parent = emptyenv())
  genv = globalenv()
  local_env[["ymat.t"]] = genv[["vanillaR.ymat.t"]]
  local_env[["ymat.c"]] = genv[["vanillaR.ymat.c"]]
  local_env[["pred.data.t"]] = genv[["vanillaR.pred.data.t"]]
  local_env[["pred.data.c"]] = genv[["vanillaR.pred.data.c"]]
  local_env[["tt"]] = genv[["vanillaR.tt"]]
  local_env[["PredictMt"]] = genv[["vanillaR.PredictMt"]]
  local_env[["PredictMc"]] = genv[["vanillaR.PredictMc"]]
  local_env[["cat.t"]] = genv[["vanillaR.cat.t"]]
  local_env[["cat.c"]] = genv[["vanillaR.cat.c"]]
  local_env[["cat.t.ctrl"]] = genv[["vanillaR.cat.t.ctrl"]]
  local_env[["cat.c.ctrl"]] = genv[["vanillaR.cat.c.ctrl"]]
  local_env[["Pr0.1"]] = genv[["vanillaR.Pr0.1"]]
  local_env[["Pr0.2"]] = genv[["vanillaR.Pr0.2"]]
  local_env[["Pr0.3"]] = genv[["vanillaR.Pr0.3"]]
  local_env[["Pr0.4"]] = genv[["vanillaR.Pr0.4"]]
  local_env[["Pr1.1"]] = genv[["vanillaR.Pr1.1"]]
  local_env[["Pr1.2"]] = genv[["vanillaR.Pr1.2"]]
  local_env[["Pr1.3"]] = genv[["vanillaR.Pr1.3"]]
  local_env[["Pr1.4"]] = genv[["vanillaR.Pr1.4"]]
  return(local_env)
}

collect_rcpp_loop_exports <- function(){
  local_env = new.env(parent = emptyenv())
  genv = globalenv()
  local_env[["ymat.t"]] = genv[["rcpp.ymat.t"]]
  local_env[["ymat.c"]] = genv[["rcpp.ymat.c"]]
  local_env[["pred.data.t"]] = genv[["rcpp.pred.data.t"]]
  local_env[["pred.data.c"]] = genv[["rcpp.pred.data.c"]]
  local_env[["tt"]] = genv[["rcpp.tt"]]
  local_env[["PredictMt"]] = genv[["rcpp.PredictMt"]]
  local_env[["PredictMc"]] = genv[["rcpp.PredictMc"]]
  local_env[["cat.t"]] = genv[["rcpp.cat.t"]]
  local_env[["cat.c"]] = genv[["rcpp.cat.c"]]
  local_env[["cat.t.ctrl"]] = genv[["rcpp.cat.t.ctrl"]]
  local_env[["cat.c.ctrl"]] = genv[["rcpp.cat.c.ctrl"]]
  local_env[["Pr0.1"]] = genv[["rcpp.Pr0.1"]]
  local_env[["Pr0.2"]] = genv[["rcpp.Pr0.2"]]
  local_env[["Pr0.3"]] = genv[["rcpp.Pr0.3"]]
  local_env[["Pr0.4"]] = genv[["rcpp.Pr0.4"]]
  local_env[["Pr1.1"]] = genv[["rcpp.Pr1.1"]]
  local_env[["Pr1.2"]] = genv[["rcpp.Pr1.2"]]
  local_env[["Pr1.3"]] = genv[["rcpp.Pr1.3"]]
  local_env[["Pr1.4"]] = genv[["rcpp.Pr1.4"]]
  return(local_env)
}

compare_values_as_vectors <- function(v1, v2){
  V1 = as.vector(as.matrix(v1))
  V2 = as.vector(as.matrix(v2))
  l1 = length(V1)
  l2 = length(V2)
  if(l1 != l2){
    warning(paste("mismatched lengths:", l1, ",", l2))
    return(F)
  }
  comp_results = vector(mode="logical", length=l1)
  for(i in (1:l1)){
    if(V1[[i]] != V2[[i]]){
      warning(paste(V1[[i]]," not equal to ",V2[[i]], sep=""))
      return(F)
    }
  }
  return(T)
}

compare_loop_envs <- function(env1, env2) {
  if(!identical(names(env1), names(env2))){
    warning("name vectors don't match between environments")
    return(F)
  }
  for(key in names(env1)){
    if(!compare_values_as_vectors(as.vector(env1[[key]]), as.vector(env2[[key]]))){
      warning(paste(key, "didn't match"))
      return(F)
    }
  }
  return(T)
}

mediation.test.compare <- function(n=100, sims=100){
  set.seed(1)
  param = reverseC::MediateDataGenerationParameters(n=n, nSim=1)
  data_item = reverseC::generateDataMatrix(param)[[1]]
  models = reverseC::assembleLinearModels(data_item)
  set.seed(1)
  pre_func_env = reverseC::test.mediate.vanilla(models@med.fit, models@out.fit, treat="X", mediator="M1", sims=sims,return_context=T, context_before=T)
  export_environment(pre_func_env)
  set.seed(1)
  cat("Running Vanilla R Mediation\n")
  tic("Vanilla R Mediation")
  vanilla_env = reverseC::test.mediate.vanilla(models@med.fit, models@out.fit, treat="X", mediator="M1", sims=sims,return_context=T, context_before=F, export_loop_vars=T)
  toc()
  vanilla_loop_env = reverseC::collect_vanillaR_loop_exports()
  set.seed(1)
  cat("Running Rcpp/Eigen Assisted Mediation with 1 Thread\n")
  tic("Rcpp/Eigen Assisted Mediation with 1 Thread")
  rcpp_env_n1 = reverseC::test.mediate.rcpp(models@med.fit, models@out.fit, treat="X", mediator="M1", sims=sims,return_context=T, context_before=F, num_cores=1, export_loop_vars=T)
  toc()
  rcpp_loop_env = reverseC::collect_rcpp_loop_exports()
  
  # rcpp_env = reverseC::test.mediate.rcpp(models@med.fit, models@out.fit, treat="X", mediator="M1", sims=100,return_context=T, context_before=F, num_cores=7)
  genv = globalenv()
  genv[["rcpp_env"]] = rcpp_env
  genv[["vanilla_env"]] = vanilla_env
  pass = T
  if(!compare_loop_envs(vanilla_loop_env, rcpp_loop_env)){
    pass = F
    warning("loop envs did not match")
  }
  cat("Running Rcpp/Eigen Assisted Mediation with 2 Threads\n")
  set.seed(1)
  tic("Rcpp/Eigen Assisted Mediation with 2 Threads")
  rcpp_env_n2 = reverseC::test.mediate.rcpp(models@med.fit, models@out.fit, treat="X", mediator="M1", sims=sims,return_context=T, context_before=F, num_cores=2)
  toc()
  cat("Running Rcpp/Eigen Assisted Mediation with 4 Threads\n")
  set.seed(1)
  tic("Rcpp/Eigen Assisted Mediation with 4 Threads")
  rcpp_env_n4 = reverseC::test.mediate.rcpp(models@med.fit, models@out.fit, treat="X", mediator="M1", sims=sims,return_context=T, context_before=F, num_cores=4)
  toc()
  cat("Running Rcpp/Eigen Assisted Mediation with 7 Threads\n")
  set.seed(1)
  tic("Rcpp/Eigen Assisted Mediation with 7 Threads")
  rcpp_env_n7 = reverseC::test.mediate.rcpp(models@med.fit, models@out.fit, treat="X", mediator="M1", sims=sims,return_context=T, context_before=F, num_cores=7)
  toc()
  #Make sure that all the rcpp thread strategies return the same et{#} matrices
  local_env = environment()
  #Make sure that the rcpp results match the vanilla R results
  for(i in 1:4){
    et_name = paste("et",i,sep="")
    for(n in c(2,4,7)){
      env_name = paste("rcpp_env_n",n,sep="")
      if(!identical(local_env[["rcpp_env_n1"]][[et_name]], local_env[[env_name]][[et_name]])){
        warning(paste("rcpp_env_n1[[",et_name ,"]] did not match ", env_name, "[[",et_name ,"]]",sep=""))
      # } else {
      #   warning(paste("rcpp_env_n1[[",et_name ,"]] matched ", env_name, "[[",et_name ,"]]",sep=""))
      }
    }
  }
  for(i in 1:4){
    et_name = paste("et",i,sep="")
    for(n in c(1,2,4,7)){
      env_name = paste("rcpp_env_n",n,sep="")
      if(!identical(local_env[["vanilla_env"]][[et_name]], local_env[[env_name]][[et_name]])){
        warning(paste("vanilla_env[[",et_name ,"]] did not match ", env_name, "[[",et_name ,"]]",sep=""))
      # } else {
      #   warning(paste("vanilla_env[[",et_name ,"]] matched ", env_name, "[[",et_name ,"]]",sep=""))
      }
    }
  }
  return(T)
}

test_all_option_combos <- function(num_jobs=2, nJobsIter=2, nSimImai=1000, n=1000){
  if(num_jobs < 2){
    stop("num_jobs must be 2 or greater. In order to test the functionality and speed effect of multiprocessing and threads, jobs=1 is already included as a default test case against your selection.")
  }
  nSim = num_jobs * nJobsIter
  cat(paste("nSim:", nSim, "; nSimImai:",nSimImai, "; n:",n, "\n", sep=""))
  for(jobs in c(1,num_jobs)){
    for(use_cpp in c(F,T)){
      for(use_multi_processing in c(T,F)){
        #do not use 1 job with multiprocessing
        bad_case_1 = jobs == 1 && use_multi_processing
        #do not use multiple jobs without either use_cpp or use_multi_processing
        bad_case_2 = (jobs > 1) && !(use_multi_processing || use_cpp)
        if(!(bad_case_1 || bad_case_2)){
          case_name = paste(paste("use_cpp:",use_cpp,"use_multi_processing:",use_multi_processing,"jobs:",jobs), sep="")
          cat(case_name)
          cat("\n")
          reverseC::reverseMAsim(use_cpp = use_cpp, use_multi_processing = use_multi_processing, num_jobs = jobs, nSim = nSim, nSimImai = nSimImai, n=n)
          cat("\n")
        }
      }
    }
  }
}