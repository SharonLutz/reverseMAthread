
#' @export
#' @title reverseDirection
#' @description placeholder_text
#' @author Sharon Lutz
#' @param n n
#' @param nSNP nSNP
#' @param MAF MAF
#' @param gamma0 gamma0
#' @param gammaX gammaX
#' @param varM varM
#' @param beta0 beta0
#' @param betaM betaM
#' @param varY varY
#' @param nSim nSim
#' @param plot.pdf plot.pdf
#' @param plot.name plot.name
#' @param alpha_level alpha_level
#' @param SEED SEED
reverseDirection <-
  function(n=1000,nSNP=10,MAF=rep(0.04,nSNP),gamma0=0,gammaX=rep(0.1,nSNP),varM=1,beta0=0,betaM=seq(from=0,to=2,length.out=4),varY=1,nSim=500,plot.pdf=T,plot.name="reverseDirection.pdf",alpha_level=0.05,SEED=1){
    
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
    
    ################################################################################
    #matrix to save MR steiger
    ################################################################################
    #save results for type 1 error rate betaM=0 and power betaM>0
    matR <- matrix(0,ncol=11,nrow=length(betaM))
    colnames(matR) <- c("cdir","cedir","pdir","pedir","sr","cdirM","cedirM","pdirM","pedirM","srM","corMY")
    
    ################################################################################
    # cycle through the simulations
    ################################################################################  
    #cycle through the simulations
    for(ii in 1:nSim){
      pcut<-50
      if(floor(ii/pcut)==ceiling(ii/pcut)){print(paste(ii,"of",nSim,"simulations"))}
      
      ################################################################################
      # cycle through values of betaM
      ################################################################################ 
      #cycle through values of betaM
      for(bM in 1:length(betaM)){
        
        ################################################################################
        # generate data
        ################################################################################                            
        # Generate X
        X <- matrix(NA,nrow=n,ncol=nSNP)
        # Create vector of SNPs
        for(jj in 1:nSNP){
          X[,jj] <- rbinom(n,2,MAF[jj])
        }
        
        # Generate M and Y
        M <- rnorm(n,gamma0 + X%*%gammaX,sqrt(varM))
        Y <- rnorm(n,beta0 + betaM[bM]*M,sqrt(varY))
        
        matR[bM,"corMY"]<-matR[bM,"corMY"]+cor(M,Y)
        
        ################################################################################
        # Get p-values from Y~X and M~X
        ################################################################################ 
        #Vector of p-values of SNP-exposure
        p_exp<-rep(0,nSNP)
        for(pe in 1:nSNP){
          p_exp[pe]<-summary(lm(M~X[,pe]))$coef[2,4]
        }
        
        #Vector of p-values of SNP-outcome
        p_out<-rep(0,nSNP)
        for(po in 1:nSNP){
          p_out[po]<-summary(lm(Y~X[,po]))$coef[2,4]
        }
        
        #Sample sizes for p_exp, p_out
        n_exp<-n
        n_out<-n
        
        #Vector of absolute correlations for SNP-exposure
        r_exp<-rep(0,nSNP)
        for(re in 1:nSNP){
          r_exp[re]<-abs(cor(X[,re],M))
        }
        
        #Vector of absolute correlations for SNP-outcome
        r_out<-rep(0,nSNP)
        for(ro in 1:nSNP){
          r_out[ro]<-abs(cor(X[,ro],Y))
        }
        
        #r_xxo Measurememt precision of exposure, default 1
        #r_yyo Measurement precision of outcome, default 1
        
        ################################################################################
        # MR Steiger correct way
        ################################################################################ 
        #A statistical test for whether the assumption that exposure causes outcome is valid
        mrs<-TwoSampleMR::mr_steiger(p_exp, p_out, n_exp, n_out, r_exp, r_out, r_xxo = 1, r_yyo = 1)
        
        #correct_causal_direction: TRUE/FALSE 
        if(mrs$correct_causal_direction==TRUE){matR[bM,"cdir"]<-matR[bM,"cdir"]+1}
        
        #steiger_test: p-value for inference of direction 
        if(mrs$steiger_test<alpha){matR[bM,"pdir"]<-matR[bM,"pdir"]+1}
        
        #correct_causal_direction_adj: TRUE/FALSE, direction of causality for given measurement error parameters 
        if(mrs$correct_causal_direction_adj==TRUE){matR[bM,"cedir"]<-matR[bM,"cedir"]+1}
        
        #steiger_test_adj: p-value for inference of direction of causality for given measurement error parameters 
        if(mrs$steiger_test_adj<alpha){matR[bM,"pedir"]<-matR[bM,"pedir"]+1}
        
        #sensitivity_ratio: Ratio of vz1/vz0. Higher means inferred direction is less susceptible to measurement error - 
        matR[bM,"sr"]<-matR[bM,"sr"]+mrs$sensitivity_ratio
        
        ################################################################################
        # MR Steiger reversed
        ################################################################################ 
        #A statistical test for whether the assumption that exposure causes outcome is valid
        mrs<-TwoSampleMR::mr_steiger(p_out, p_exp, n_out, n_exp, r_out, r_exp, r_xxo = 1, r_yyo = 1)
        
        #correct_causal_direction: TRUE/FALSE 
        if(mrs$correct_causal_direction==TRUE){matR[bM,"cdirM"]<-matR[bM,"cdirM"]+1}
        
        #steiger_test: p-value for inference of direction 
        if(mrs$steiger_test<alpha){matR[bM,"pdirM"]<-matR[bM,"pdirM"]+1}
        
        #correct_causal_direction_adj: TRUE/FALSE, direction of causality for given measurement error parameters 
        if(mrs$correct_causal_direction_adj==TRUE){matR[bM,"cedirM"]<-matR[bM,"cedirM"]+1}
        
        #steiger_test_adj: p-value for inference of direction of causality for given measurement error parameters 
        if(mrs$steiger_test_adj<alpha){matR[bM,"pedirM"]<-matR[bM,"pedirM"]+1}
        
        #sensitivity_ratio: Ratio of vz1/vz0. Higher means inferred direction is less susceptible to measurement error - 
        matR[bM,"srM"]<-matR[bM,"srM"]+mrs$sensitivity_ratio
        
        ################################################################################
        # end loops
        ################################################################################    
      }#beta loop
    }#sim loop
    
    mat_total <- matR/nSim
    
    if(plot.pdf){
      pdf(plot.name)
      plot(-2,-2,xlim=c(min(betaM),max(betaM)+0.05),ylim=c(0.5,1),main="",xlab=expression(beta),ylab="Correct Direction (% of simulations)")
      lines(betaM,mat_total[,"cdir"],col=1,pch=1,type="b")
      dev.off()
    }
    
    # Print out the matrix but use list
    list(mat_total)
    
  }
