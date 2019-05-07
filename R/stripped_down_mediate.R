#' @include pval_func.R
#our calls always this style: mediate(model.m, model.y, treat= "X", mediator="M", sims = nSimImai)
#both models are lm with no special types / conditions
stripped_down_mediate <- function(model.m, model.y, sims = 1000, treat = "treat.name", mediator = "med.name",
                                  conf.level = .95, control.value = 0, treat.value = 1){
  
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
  
  effects.tmp <- array(NA, dim = c(n, sims, 4))
  
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
      
      Pr1[,j] <- t(as.matrix(YModel[j,])) %*% t(ymat.t)
      Pr0[,j] <- t(as.matrix(YModel[j,])) %*% t(ymat.c)
      
      rm(ymat.t, ymat.c, pred.data.t, pred.data.c)
    }
    
    effects.tmp[,,e] <- Pr1 - Pr0 ### e=1:mediation(1); e=2:mediation(0); e=3:direct(1); e=4:direct(0)
    rm(Pr1, Pr0)
  }
  
  rm(PredictM1, PredictM0, YModel, MModel)
  
  et1<-effects.tmp[,,1] ### mediation effect (1)
  et2<-effects.tmp[,,2] ### mediation effect (0)
  et3<-effects.tmp[,,3] ### direct effect (1)
  et4<-effects.tmp[,,4] ### direct effect (0)
  
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
  d0.p <- pval(delta.0, d0)
  d1.p <- pval(delta.1, d1)
  d.avg.p <- pval(delta.avg, d.avg)
  z0.p <- pval(zeta.0, z0)
  z1.p <- pval(zeta.1, z1)
  z.avg.p <- pval(zeta.avg, z.avg)        
  n0.p <- pval(nu.0, n0)
  n1.p <- pval(nu.1, n1)
  n.avg.p <- pval(nu.avg, n.avg)
  tau.p <- pval(tau, tau.coef)
  
  # Detect whether models include T-M interaction
  INT <- paste(treat,mediator,sep=":") %in% attr(terms(model.y),"term.labels") |
    paste(mediator,treat,sep=":") %in% attr(terms(model.y),"term.labels")
  
  return(SimpleMediateResult(direct_p = z.avg.p, indirect_p=d.avg.p))
}
