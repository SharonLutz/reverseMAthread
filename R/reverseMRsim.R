#' @importFrom MendelianRandomization mr_input mr_egger mr_ivw mr_median
 
#' @export
reverseMRsim <-
function(n=1000,nSNP=3,MAF=c(0.2,0.2,0.2),gamma0=0,gammaX=c(0.2,0.2,0.2),varM=1,beta0=0,betaM=c(0,0.1),
                         varY=1,nSim=100,plot.pdf=T,plot.name="reverseMRsim.pdf",alpha_level=0.05,SEED=1){
  # Set the seed
  set.seed(SEED)
  
  # Error checks
  if(n<0 | n==0 | floor(n)!=ceiling(n) ){stop("Error: n must be an integer greater than or equal to 1")}
  if(nSNP<0 | nSNP==0 | floor(nSNP)!=ceiling(nSNP) ){stop("Error: nSNP must be an integer greater than or equal to 1")}
  if(length(MAF)!=nSNP){stop("Error: nSNP must equal the length of MAF")}
  if(min(MAF)<0 | max(MAF)>1){stop("Error: MAF must be greater than 0 and less than 1")}
  if(length(unique(betaM))<2){stop("Error: betaM must be a vector with at least two values")}
  if(nSNP<3){stop("Error: nSNP must be at least 3")}
  if(!varM>0){stop("Error: varM must be greater than 0")}
  if(!varY>0){stop("Error: varY must be greater than 0")}
  if(length(gammaX)!=nSNP){stop("Error: ")}
  
  # Create total matrix
  mat_total <- matrix(0,nrow=6,ncol=length(betaM))
  rownames(mat_total) <- c("PMR.EggerNR","PMR.IVWNR","PMR.MedianNR","PMR.EggerR","PMR.IVWR","PMR.MedianR")

  for(i in 1:nSim){
    if(floor(i/5)==ceiling(i/5)){print(paste(i,"of",nSim,"simulations"))}
    
    # Creat a results matrix for each simulation
    mat_results <- matrix(0,nrow=6,ncol=length(betaM))
    rownames(mat_results) <- c("PMR.EggerNR","PMR.IVWNR","PMR.MedianNR","PMR.EggerR","PMR.IVWR","PMR.MedianR")
    
    for(bM.ind in 1:length(betaM)){
      
      # Generate X
      X <- matrix(NA,nrow=n,ncol=nSNP)
      # Create vector of SNPs
      for(jj in 1:nSNP){
        X[,jj] <- rbinom(n,2,MAF[jj])
      }
      
      # Generate M and Y
      M <- rnorm(n,gamma0 + X%*%gammaX,sqrt(varM))
      Y <- rnorm(n,beta0 + betaM[bM.ind]*M,sqrt(varY))
      
      # Get the input for MR methods
      betaXM <- c(); betaXMse <- c()
      betaXY <- c(); betaXYse <- c()
      
      for(snp in 1:nSNP){
        # betaX are the beta-coefficients from univariable regression analyses of the exposure on each
        # genetic variant in turn, and betaXse are the standard errors
        fitM <- lm(M~X[,snp])
        betaXM <- c(betaXM,summary(fitM)$coefficients[2,1])
        betaXMse <- c(betaXMse,summary(fitM)$coefficients[2,2])
        
        # betaXy are the beta-coefficients from regression analyses of the outcome on each genetic
        # variant in turn, and betaXsey are the standard errors
        fitY <- lm(Y~X[,snp])
        betaXY <- c(betaXY,summary(fitY)$coefficients[2,1])
        betaXYse <- c(betaXYse,summary(fitY)$coefficients[2,2])
      }
      
      # Create MR.input object
      mr.input <- mr_input(bx = betaXM, bxse = betaXMse, by = betaXY, byse = betaXYse)
      mr.input2 <- mr_input(bx = betaXY, bxse = betaXYse, by = betaXM, byse = betaXMse)
      
      # Get the estimates and p-values from the methods
      if(mr_egger(mr.input,robust = F,penalized = F)$Pvalue.Est < alpha_level){
        mat_results["PMR.EggerNR",bM.ind] <- mat_results["PMR.EggerNR",bM.ind] +1
      }
      if(mr_ivw(mr.input,robust = F,penalized = F)$Pvalue < alpha_level){
        mat_results["PMR.IVWNR",bM.ind] <- mat_results["PMR.IVWNR",bM.ind] +1
      }
      if(mr_median(mr.input,seed = NA)$Pvalue < alpha_level){
        mat_results["PMR.MedianNR",bM.ind] <- mat_results["PMR.MedianNR",bM.ind] +1
      }
      
      # Get the estimates and p-values from the methods reversed
      if(mr_egger(mr.input2,robust = F,penalized = F)$Pvalue.Est < alpha_level){
        mat_results["PMR.EggerR",bM.ind] <- mat_results["PMR.EggerR",bM.ind] +1
      }
      if(mr_ivw(mr.input2,robust = F,penalized = F)$Pvalue < alpha_level){
        mat_results["PMR.IVWR",bM.ind] <- mat_results["PMR.IVWR",bM.ind] +1
      }
      if(mr_median(mr.input2,seed = NA)$Pvalue < alpha_level){
        mat_results["PMR.MedianR",bM.ind] <- mat_results["PMR.MedianR",bM.ind] +1
      }
      
    } # End of betaM loop
    
    mat_total <- mat_total + mat_results
    
  } # End of nSim
  
  mat_total <- mat_total/nSim
  
  if(plot.pdf){
    pdf(plot.name)
    plot(-1,-1, xlim=c(0,max(betaM)), ylim=c(0,1),xlab="betaM values",ylab="")
    points(betaM,mat_total["PMR.EggerNR",],type="b",lty=2,col=1,pch=1)
    points(betaM,mat_total["PMR.IVWNR",],type="b",lty=3,col=2,pch=2)
    points(betaM,mat_total["PMR.MedianNR",],type="b",lty=4,col=3,pch=3)
    points(betaM,mat_total["PMR.EggerR",],type="b",lty=5,col=4,pch=4)
    points(betaM,mat_total["PMR.IVWR",],type="b",lty=6,col=5,pch=5)
    points(betaM,mat_total["PMR.MedianR",],type="b",lty=2,col=6,pch=6)
    legend("topleft",lty=c(2:6),col=c(1:6),pch=c(1:6),
           legend=c("PMR.EggerNR","PMR.IVWNR","PMR.MedianNR","PMR.EggerR","PMR.IVWR","PMR.MedianR"),cex = 0.6)
    dev.off()
  }
 
  # Print out the matrix but use list
  list(mat_total)
  
}
