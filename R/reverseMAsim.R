reverseMAsim <-
function(n=1000,pX=0.2,gamma0=0,gammaX=0.1,varM=1,beta0=0,betaX=1,betaM=c(0,0.1,0.2),varY=1,
                         nSim=100,nSimImai=1000,SEED=1,plot.pdf=T,plot.name="reverseMAsim.pdf",alpha_level=0.05){
  library(mediation)
  
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
  
  for(i in 1:nSim){
    if(floor(i/10)==ceiling(i/10)){print(paste(i,"of",nSim,"simulations"))}
      
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
      med.out <- mediate(med.fit,out.fit,treat = "X",mediator = "M",sims = nSimImai)
      
      # Get the direct and indirect effects
      pval_direct <- summary(med.out)$z.avg.p
      pval_indirect <- summary(med.out)$d.avg.p
      
      # Add to the matrix
      if(pval_direct<alpha_level){mat_results[bM.ind,"DirectNR"] <- mat_results[bM.ind,"DirectNR"]+1 }
      if(pval_indirect<alpha_level){mat_results[bM.ind,"IndirectNR"] <- mat_results[bM.ind,"IndirectNR"]+1 }
      
      # Reverse the order and run the mediation model again
      M2 <- Y
      Y2 <- M
      
      # Fit the mediation model
      med.fitR <- (lm(M2~X))
      out.fitR <- (lm(Y2~X+M2))
      med.outR <- mediate(med.fitR,out.fitR,treat = "X",mediator = "M2",sims = nSimImai)
      
      # Get the direct and indirect effects
      pval_direct_r <- summary(med.outR)$z.avg.p
      pval_indirect_r <- summary(med.outR)$d.avg.p
      
      # Add to the matrix
      if(pval_direct_r<alpha_level){mat_results[bM.ind,"DirectR"] <- mat_results[bM.ind,"DirectR"]+1 }
      if(pval_indirect_r<alpha_level){mat_results[bM.ind,"IndirectR"] <- mat_results[bM.ind,"IndirectR"]+1 }
      
    } # End of bM.ind
    
    mat_total <- mat_total+mat_results
    
  } # End of nSim
  
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
