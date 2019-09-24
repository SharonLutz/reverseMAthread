#' @include pval.R
#our calls always this style: mediate(model.m, model.y, treat= "X", mediator="M", sims = nSimImai)
#both models are lm with no special types / conditions

mediate_with_rcpp <- 
  function(model.m, model.y, sims = 1000,boot = FALSE, boot.ci.type = "perc", treat = "treat.name", mediator = "med.name",covariates = NULL,
           conf.level = .95, control.value = 0, treat.value = 1, long = TRUE,robustSE = FALSE, cluster = NULL){
    
    cl <- match.call()
    
    num_threads = getOption("mediate.threads", default = 1)
    
    # Model type indicators
    isGam.y <- inherits(model.y, "gam")
    isGam.m <- inherits(model.m, "gam")
    isGlm.y <- inherits(model.y, "glm")  # Note gam and bayesglm also inherits "glm"
    isGlm.m <- inherits(model.m, "glm")  # Note gam and bayesglm also inherits "glm"
    isLm.y <- inherits(model.y, "lm")    # Note gam, glm and bayesglm also inherit "lm"
    isLm.m <- inherits(model.m, "lm")    # Note gam, glm and bayesglm also inherit "lm"
    isVglm.y <- inherits(model.y, "vglm")
    isRq.y <- inherits(model.y, "rq")
    isRq.m <- inherits(model.m, "rq")
    isOrdered.y <- inherits(model.y, "polr")  # Note bayespolr also inherits "polr"
    isOrdered.m <- inherits(model.m, "polr")  # Note bayespolr also inherits "polr"
    isSurvreg.y <- inherits(model.y, "survreg")
    isSurvreg.m <- inherits(model.m, "survreg")
    isMer.y <- inherits(model.y, "merMod") # Note lmer and glmer do not inherit "lm" and "glm"
    isMer.m <- inherits(model.m, "merMod") # Note lmer and glmer do not inherit "lm" and "glm"
    
    if(isGlm.m){
      FamilyM <- model.m$family$family
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
    
    if(!is.null(weights.m) && isGlm.m && FamilyM == "binomial"){
      message("weights taken as sampling weights, not total number of trials")
    }
    
    if(is.null(weights.m)){
      weights.m <- rep(1,nrow(m.data))
    }
    if(is.null(weights.y)){
      weights.y <- rep(1,nrow(y.data))
    }
    
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
    
    if(isGlm.m){
      muM1 <- model.m$family$linkinv(tcrossprod(MModel, mmat.t))
      muM0 <- model.m$family$linkinv(tcrossprod(MModel, mmat.c))
      
      if(FamilyM == "poisson"){
        PredictM1 <- matrix(rpois(sims*n, lambda = muM1), nrow = sims)
        PredictM0 <- matrix(rpois(sims*n, lambda = muM0), nrow = sims)
      } else if (FamilyM == "Gamma") {
        shape <- gamma.shape(model.m)$alpha
        PredictM1 <- matrix(rgamma(n*sims, shape = shape,
                                   scale = muM1/shape), nrow = sims)
        PredictM0 <- matrix(rgamma(n*sims, shape = shape,
                                   scale = muM0/shape), nrow = sims)
      } else if (FamilyM == "binomial"){
        PredictM1 <- matrix(rbinom(n*sims, size = 1,
                                   prob = muM1), nrow = sims)
        PredictM0 <- matrix(rbinom(n*sims, size = 1,
                                   prob = muM0), nrow = sims)
      } else if (FamilyM == "gaussian"){
        sigma <- sqrt(summary(model.m)$dispersion)
        error <- rnorm(sims*n, mean=0, sd=sigma)
        PredictM1 <- muM1 + matrix(error, nrow=sims)
        PredictM0 <- muM0 + matrix(error, nrow=sims)
      } else if (FamilyM == "inverse.gaussian"){
        disp <- summary(model.m)$dispersion
        PredictM1 <- matrix(SuppDists::rinvGauss(n*sims, nu = muM1,
                                                 lambda = 1/disp), nrow = sims)
        PredictM0 <- matrix(SuppDists::rinvGauss(n*sims, nu = muM0,
                                                 lambda = 1/disp), nrow = sims)
      } else {
        stop("unsupported glm family")
      }
      
      ### Case I-1-c: Linear
    } else if(isLm.m){
      sigma <- summary(model.m)$sigma
      error <- rnorm(sims*n, mean=0, sd=sigma)
      muM1 <- tcrossprod(MModel, mmat.t)
      muM0 <- tcrossprod(MModel, mmat.c)
      PredictM1 <- muM1 + matrix(error, nrow=sims)
      PredictM0 <- muM0 + matrix(error, nrow=sims)
      rm(error)
    }
    
    rm(mmat.t, mmat.c)
    
    #####################################
    ##  Outcome Predictions
    #####################################
    
    if(num_threads > 1){
      threaded_mediate_helper(environment(), num_threads)
    } else {
      mediate_helper(environment())
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
    d0.p <- mediate_pval(delta.0, d0)
    d1.p <- mediate_pval(delta.1, d1)
    d.avg.p <- mediate_pval(delta.avg, d.avg)
    z0.p <- mediate_pval(zeta.0, z0)
    z1.p <- mediate_pval(zeta.1, z1)
    z.avg.p <- mediate_pval(zeta.avg, z.avg)        
    n0.p <- mediate_pval(nu.0, n0)
    n1.p <- mediate_pval(nu.1, n1)
    n.avg.p <- mediate_pval(nu.avg, n.avg)
    tau.p <- mediate_pval(tau, tau.coef)
    
    # Detect whether models include T-M interaction
    INT <- paste(treat,mediator,sep=":") %in% attr(terms(model.y),"term.labels") |
      paste(mediator,treat,sep=":") %in% attr(terms(model.y),"term.labels")
    
    if(long && !isMer.y && !isMer.m) {
      out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                  d0.p=d0.p, d1.p=d1.p,
                  d0.sims=delta.0, d1.sims=delta.1,
                  z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
                  z0.p=z0.p, z1.p=z1.p,
                  z0.sims=zeta.0, z1.sims=zeta.1,
                  n0=n0, n1=n1, n0.ci=n0.ci, n1.ci=n1.ci,
                  n0.p=n0.p, n1.p=n1.p,
                  n0.sims=nu.0, n1.sims=nu.1,
                  tau.coef=tau.coef, tau.ci=tau.ci, tau.p=tau.p,
                  tau.sims=tau,
                  d.avg=d.avg, d.avg.p=d.avg.p, d.avg.ci=d.avg.ci, d.avg.sims=delta.avg,
                  z.avg=z.avg, z.avg.p=z.avg.p, z.avg.ci=z.avg.ci, z.avg.sims=zeta.avg,
                  n.avg=n.avg, n.avg.p=n.avg.p, n.avg.ci=n.avg.ci, n.avg.sims=nu.avg,
                  boot=boot, boot.ci.type=boot.ci.type,
                  treat=treat, mediator=mediator,
                  covariates=covariates,
                  INT=INT, conf.level=conf.level,
                  model.y=model.y, model.m=model.m,
                  control.value=control.value, treat.value=treat.value,
                  nobs=n, sims=sims, call=cl,
                  robustSE = robustSE, cluster = cluster)
      class(out) <- "mediate"
    } 
    if(!long && !isMer.y && !isMer.m){
      out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                  d0.p=d0.p, d1.p=d1.p,
                  z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
                  z0.p=z0.p, z1.p=z1.p,
                  n0=n0, n1=n1, n0.ci=n0.ci, n1.ci=n1.ci,
                  n0.p=n0.p, n1.p=n1.p,
                  tau.coef=tau.coef, tau.ci=tau.ci, tau.p=tau.p,
                  d.avg=d.avg, d.avg.p=d.avg.p, d.avg.ci=d.avg.ci,
                  z.avg=z.avg, z.avg.p=z.avg.p, z.avg.ci=z.avg.ci,
                  n.avg=n.avg, n.avg.p=n.avg.p, n.avg.ci=n.avg.ci,
                  boot=boot, boot.ci.type=boot.ci.type,
                  treat=treat, mediator=mediator,
                  covariates=covariates,
                  INT=INT, conf.level=conf.level,
                  model.y=model.y, model.m=model.m,
                  control.value=control.value, treat.value=treat.value,
                  nobs=n, sims=sims, call=cl,
                  robustSE = robustSE, cluster = cluster)
      class(out) <- "mediate"
    }
    return(out)
  }