#' @include mediate_s4_classes.R multi_process_mediate.R

#' @title reverseMAsim.SingleProcess
#' @description A function to similate the performance of the mediate function from the mediation package in scenarios of reverse causality. Leverages threading and Rcpp with Eigen for faster computation speed.
#' @author Michael Gooch, Annie Thwing, Sharon Lutz
#' @param n is the sample size.
#' @param pX is the minor allele frequency
#' @param gamma0 is the intercept for M
#' @param gammaX is the association of X with M
#' @param varM is the variance of M
#' @param beta0 is the intercept for Y
#' @param betaX is the direct effect of X on Y
#' @param betaM is a vector of different associations of M with Y
#' @param varY is the variance of Y
#' @param nSim is the number of simulations to run
#' @param nSimImai is the number of simulations to run in mediate from the mediation package
#' @param SEED is the seed
#' @param plot.pdf is T to output a plot, is F to not output a plot
#' @param plot.name is the name of the plot
#' @param alpha_level is the significance level
#' @param use_cpp use Rcpp with Eigen
#' @param num_jobs the number of cores to use, i.e. the number of threads to spawn
#' @return a matrix of the power of the mediate method from the mediation package to detect an effect of the mediator M on the outcome Y when M and Y are correctly specified and also when they are incorrectly specified (the true mediator is Y and the true outcome is M)
reverseMAsim.SingleProcess <- 
  function(n=1000,pX=0.2,gamma0=0,gammaX=0.1,varM=1,beta0=0,betaX=1,betaM=c(0,0.1,0.2),varY=1,
           nSim=100,nSimImai=1000,SEED=1,plot.pdf=T,plot.name="reverseMAsim.pdf",alpha_level=0.05, 
           use_cpp=F, num_jobs=1){
    
    # cat("using One Process\n")
    
    if((!use_cpp) && (num_jobs > 1)){
      stop(paste("num_jobs=",num_jobs, " will have no effect when use_cpp=F and use_multiprocessing=F",sep=""))
    }
    options(mediate.jobs = 1)
    options(mediate.threads = num_jobs)
    
    # Set the seed.
    set.seed(SEED)
    
    next_seed = SEED + 1
    g_env = globalenv()
    data_gen_randstate = g_env[[".Random.seed"]]
    # Error checks
    if(n<=0 || floor(n)!=ceiling(n) ){stop("Error: n must be an integer greater than or equal to 1")}
    if(pX<0 | pX>1){stop("Error: pX must be greater than 0 and less than 1")}
    if(varM<=0){stop("Error: varM must be greater than 0")}
    if(varY<=0){stop("Error: varY must be greater than 0")}
    if(length(unique(betaM))!=length(betaM)){stop("Error: betaM must contain unique values")}
    if(length(betaM)<2){stop("Error: betaM must be a vector with at least two values")}
    if(alpha_level>1 | alpha_level<0){stop("Error: alpha_level must be between 0 and 1")}
    
    mat_total <- matrix(0,nrow=length(betaM),ncol=4)
    colnames(mat_total) <- c("DirectNR","IndirectNR","DirectR","IndirectR")
    pb_opts = getOption("pboptions")
    pb = pbapply::timerProgressBar(min=0, max=nSim, style=pb_opts$style, char = pb_opts$char, width = pb_opts$txt.width, title = pb_opts$title, label=pb_opts$label)
    for(i in 1:nSim){
      
      # Create matrix to store the results
      mat_results <- matrix(0,nrow=length(betaM),ncol=4)
      colnames(mat_results) <- c("DirectNR","IndirectNR","DirectR","IndirectR")
      
      for(bM.ind in 1:length(betaM)){
        # Generate the data
        
        g_env[[".Random.seed"]] = data_gen_randstate
        
        X <- rbinom(n,2,pX)
        M <- rnorm(n,gamma0 + gammaX*X,sqrt(varM))
        Y <- rnorm(n,beta0 + betaX*X + betaM[bM.ind]*M,sqrt(varY))
        
        data_gen_randstate = g_env[[".Random.seed"]]
        
        set.seed(next_seed)
        next_seed = next_seed + 1
        # Fit the mediation model
        med.fit <- (lm(M~X))
        out.fit <- (lm(Y~X+M))
        
        # Reverse the order and run the mediation model again
        M2 <- Y
        Y2 <- M
        
        # Fit the mediation model
        med.fitR <- (lm(M2~X))
        out.fitR <- (lm(Y2~X+M2))
        
        pval_direct=0.0
        pval_indirect=0.0
        pval_direct_r=0.0
        pval_indirect_r=0.0
        
        if(!use_cpp){
          med.out <- mediation::mediate(med.fit,out.fit,treat = "X",mediator = "M",sims = nSimImai)
          med.outR <- mediation::mediate(med.fitR,out.fitR,treat = "X",mediator = "M2",sims = nSimImai)
          # Get the direct and indirect effects
          pval_direct <- summary(med.out)$z.avg.p
          pval_indirect <- summary(med.out)$d.avg.p
          # Get the direct and indirect effects
          pval_direct_r <- summary(med.outR)$z.avg.p
          pval_indirect_r <- summary(med.outR)$d.avg.p
        } else {
          med.out <- reverseC::mediate_with_rcpp (med.fit,out.fit,treat = "X",mediator = "M",sims = nSimImai)
          med.outR <- reverseC::mediate_with_rcpp(med.fitR,out.fitR,treat = "X",mediator = "M2",sims = nSimImai)
          # Get the direct and indirect effects
          pval_direct <- med.out@direct_p
          pval_indirect <- med.out@indirect_p
          # Get the direct and indirect effects
          pval_direct_r <- med.outR@direct_p
          pval_indirect_r <- med.outR@indirect_p
        }
        
        
        # Add to the matrix
        if(pval_direct<alpha_level){mat_results[bM.ind,"DirectNR"] <- mat_results[bM.ind,"DirectNR"]+1 }
        if(pval_indirect<alpha_level){mat_results[bM.ind,"IndirectNR"] <- mat_results[bM.ind,"IndirectNR"]+1 }
        
        # Add to the matrix
        if(pval_direct_r<alpha_level){mat_results[bM.ind,"DirectR"] <- mat_results[bM.ind,"DirectR"]+1 }
        if(pval_indirect_r<alpha_level){mat_results[bM.ind,"IndirectR"] <- mat_results[bM.ind,"IndirectR"]+1 }
        
      } # End of bM.ind
      
      mat_total <- mat_total+mat_results
      pbapply::setTimerProgressBar(pb, value=i)
    } # End of nSim
    
    pbapply::closepb(pb)
    mat_total <- mat_total/nSim
    
    if(plot.pdf){
      pdf(plot.name)
      plot(-1,-1, xlim=c(min(betaM),max(betaM)), ylim=c(0,1),xlab="betaM values",ylab="")
      points(betaM,mat_total[,"DirectNR"],type="b",lty=2,col=1,pch=1)
      points(betaM,mat_total[,"IndirectNR"],type="b",lty=3,col=2,pch=2)
      points(betaM,mat_total[,"DirectR"],type="b",lty=4,col=3,pch=3)
      points(betaM,mat_total[,"IndirectR"],type="b",lty=5,col=4,pch=4)
      legend("left",lty=c(2:5),col=c(1:4),pch=c(1:4),legend=c("DirectNR","IndirectNR","DirectR","IndirectR"))
      dev.off()
    }
    
    # Print out the matrix
    list(mat_total)
  }

#' @title reverseMAsim.MultiProcess
#' @description A function to similate the performance of the mediate function from the mediation package in scenarios of reverse causality. Leverages parallel processing and Rcpp with Eigen for faster computation speed.
#' @author Michael Gooch, Annie Thwing, Sharon Lutz
#' @param n is the sample size.
#' @param pX is the minor allele frequency
#' @param gamma0 is the intercept for M
#' @param gammaX is the association of X with M
#' @param varM is the variance of M
#' @param beta0 is the intercept for Y
#' @param betaX is the direct effect of X on Y
#' @param betaM is a vector of different associations of M with Y
#' @param varY is the variance of Y
#' @param nSim is the number of simulations to run
#' @param nSimImai is the number of simulations to run in mediate from the mediation package
#' @param SEED is the seed
#' @param plot.pdf is T to output a plot, is F to not output a plot
#' @param plot.name is the name of the plot
#' @param alpha_level is the significance level
#' @param use_cpp use Rcpp with Eigen
#' @param num_jobs the number of cores to use, i.e. the number of parallel procesess to spawn
#' @return a matrix of the power of the mediate method from the mediation package to detect an effect of the mediator M on the outcome Y when M and Y are correctly specified and also when they are incorrectly specified (the true mediator is Y and the true outcome is M)
reverseMAsim.MultiProcess <- 
  function(n=1000,pX=0.2,gamma0=0,gammaX=0.1,varM=1,beta0=0,betaX=1,betaM=c(0,0.1,0.2),varY=1,
           nSim=100,nSimImai=1000,SEED=1,plot.pdf=T,plot.name="reverseMAsim.pdf",alpha_level=0.05, 
           use_cpp=F, num_jobs=2){
    # cat("using MultiProcessing\n")
    # Set the seed.
    set.seed(SEED)
    
    if(parallel::detectCores() == 1){
      warning("your machine may not be suitable for multiprocessing, only 1 core was detected")
    }
    if(num_jobs < 2){
      stop("There is no point in using MultiProcessing with less than 2 jobs")
    }
    
    if((nSim / num_jobs) < 1.0){
      warning(paste("you don't have enough Simulations in nSim:",nSim," to fully benefit from num_jobs:",num_jobs,sep=""))
    }
    options(mediate.jobs = num_jobs)
    options(mediate.threads = 1)
    
    if(alpha_level>1 | alpha_level<0){stop("Error: alpha_level must be between 0 and 1")}
    if(length(nSimImai) != 1){stop ("Error: nSimImai must be a single integer value")}
    if(nSimImai<=0 || floor(nSimImai)!=ceiling(nSimImai) ){stop("Error: n must be an integer greater than or equal to 1")}
    
    #other values will be checked by the validator of the class
    med_vars = MediateDataGenerationParameters(
      n = n,
      pX = pX,
      gamma0 = gamma0,
      gammaX = gammaX,
      varM = varM,
      beta0 = beta0,
      betaX = betaX,
      betaM = betaM,
      varY = varY,
      nSim = nSim
    )
    
    mat_total <- matrix(0,nrow=length(betaM),ncol=4)
    colnames(mat_total) <- c("DirectNR","IndirectNR","DirectR","IndirectR")
    #generate the data needed to make linear models.
    data_matrix = generateDataMatrix(med_vars, SEED)
    # cat("running mediation on models")
    
    result.matrix = mediate_parallel(data_matrix, nSimImai)
    
    rm(data_matrix)
    
    # print("processing results")
    for(i in 1:nSim){
      #if(floor(i/10)==ceiling(i/10)){print(paste(i,"of",nSim,"simulations"))}
      
      # Create matrix to store the results
      mat_results <- matrix(0,nrow=length(betaM),ncol=4)
      colnames(mat_results) <- c("DirectNR","IndirectNR","DirectR","IndirectR")
      
      for(bM.ind in 1:length(betaM)){
        
        # Get the direct and indirect effects
        pval_direct <- result.matrix[[bM.ind,i]]@pval_direct
        pval_indirect <- result.matrix[[bM.ind,i]]@pval_indirect
        
        # Add to the matrix
        if(pval_direct<alpha_level){mat_results[bM.ind,"DirectNR"] <- mat_results[bM.ind,"DirectNR"]+1 }
        if(pval_indirect<alpha_level){mat_results[bM.ind,"IndirectNR"] <- mat_results[bM.ind,"IndirectNR"]+1 }
        
        # Get the direct and indirect effects
        pval_direct_r <- result.matrix[[bM.ind,i]]@pval_direct_r
        pval_indirect_r <- result.matrix[[bM.ind,i]]@pval_indirect_r
        
        # Add to the matrix
        if(pval_direct_r<alpha_level){mat_results[bM.ind,"DirectR"] <- mat_results[bM.ind,"DirectR"]+1 }
        if(pval_indirect_r<alpha_level){mat_results[bM.ind,"IndirectR"] <- mat_results[bM.ind,"IndirectR"]+1 }
        
      } # End of bM.ind
      
      mat_total <- mat_total+mat_results
      
    } # End of nSim
    
    rm(result.matrix)
    
    mat_total <- mat_total/nSim
    
    if(plot.pdf){
      pdf(plot.name)
      plot(-1,-1, xlim=c(min(betaM),max(betaM)), ylim=c(0,1),xlab="betaM values",ylab="")
      points(betaM,mat_total[,"DirectNR"],type="b",lty=2,col=1,pch=1)
      points(betaM,mat_total[,"IndirectNR"],type="b",lty=3,col=2,pch=2)
      points(betaM,mat_total[,"DirectR"],type="b",lty=4,col=3,pch=3)
      points(betaM,mat_total[,"IndirectR"],type="b",lty=5,col=4,pch=4)
      legend("left",lty=c(2:5),col=c(1:4),pch=c(1:4),legend=c("DirectNR","IndirectNR","DirectR","IndirectR"))
      dev.off()
    }
    
    # Print out the matrix
    list(mat_total)
  }
