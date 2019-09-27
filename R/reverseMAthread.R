#' @include generate_data_matrix.R
#' @include input_error_check.R
#' @include simulate_and_mediate.R
#' @include multi_process_mediate.R
#' @import mediation
#' @import pbapply

#' @export
reverseMAthread <-
  function(n=1000,pX=0.2,gamma0=0,gammaX=0.1,varM=1,beta0=0,betaX=1,betaM=c(0,0.1,0.2),varY=1,
           nSim=100,nSimImai=1000,SEED=1,plot.pdf=T,plot.name="reverseMAsim.pdf",alpha_level=0.05, 
           use_multi_processing=F, use_cpp=F, num_jobs=1){
    
    pbapply::pboptions(style=1,type="timer")
    
    input_error_check(n=n,pX=pX,gamma0=gamma0,gammaX=gammaX,varM=varM,beta0=beta0,betaX=betaX,betaM=betaM,varY=varY,
                      nSim=nSim,nSimImai=nSimImai,alpha_level=alpha_level, 
                      use_multi_processing=use_multi_processing, use_cpp=use_cpp, num_jobs=num_jobs)
    
    mat_total <- matrix(0,nrow=length(betaM),ncol=4)
    colnames(mat_total) <- c("DirectNR","IndirectNR","DirectR","IndirectR")
    #generate the data needed to make linear models.
    data_matrix = generate_data_matrix(n=n, pX=pX, gamma0=gamma0, gammaX=gammaX, varM=varM, 
                                       beta0=beta0, betaX=betaX, betaM=betaM, varY=varY, nSim=nSim, 
                                       nSimImai=nSimImai, SEED=SEED, use_cpp=use_cpp, num_jobs=num_jobs)
    
    # cat("running mediation on models")
    
    if(use_multi_processing){
      result.matrix = mediate_parallel(data_matrix, num_jobs = num_jobs)
    } else {
      result.matrix = pbapply::pblapply(data_matrix, simulate_and_mediate)
      dim(result.matrix) = dim(data_matrix)
    }
    rm(data_matrix)
    
    # print("processing results")
    for(i in 1:nSim){
      #if(floor(i/10)==ceiling(i/10)){print(paste(i,"of",nSim,"simulations"))}
      
      # Create matrix to store the results
      mat_results <- matrix(0,nrow=length(betaM),ncol=4)
      colnames(mat_results) <- c("DirectNR","IndirectNR","DirectR","IndirectR")
      
      for(bM.ind in 1:length(betaM)){
        data_matrix_element = result.matrix[[bM.ind,i]]
        
        summary_obj = mediation::summary.mediate(data_matrix_element[["med.out"]])
        summary_obj.r = mediation::summary.mediate(data_matrix_element[["med.out.r"]])
        
        # Get the direct and indirect effects
        pval_direct <- summary_obj[["z.avg.p"]]
        pval_indirect <- summary_obj[["d.avg.p"]]
        
        # Add to the matrix
        if(pval_direct<alpha_level){mat_results[bM.ind,"DirectNR"] <- mat_results[bM.ind,"DirectNR"]+1 }
        if(pval_indirect<alpha_level){mat_results[bM.ind,"IndirectNR"] <- mat_results[bM.ind,"IndirectNR"]+1 }
        
        # Get the direct and indirect effects
        pval_direct_r <- summary_obj.r[["z.avg.p"]]
        pval_indirect_r <- summary_obj.r[["d.avg.p"]]
        
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