#' @importFrom mediation mediate

#' @export
reverseMAdata <-
function(X, M, Y, SEED=1,nSimImai=1000){
  
  # Set the seed.
  set.seed(SEED)
  
  # Error checks
  if(length(X)!=length(M)){stop("Error: the rows of X must equal the length of M")}
  if(length(X)!=length(Y)){stop("Error: the rows of X must equal the length of Y")}
  # if(nrow(X)!=length(M)){stop("Error: the rows of X must equal the length of M")}
  # if(nrow(X)!=length(Y)){stop("Error: the rows of X must equal the length of Y")}
  
  # Create matrix to store the results
  mat_results <- matrix(NA,nrow=8,ncol=1)
  rownames(mat_results) <- c("EDirectNR","PDirectNR","EIndirectNR","PIndirectNR","EDirectR","PDirectR","EIndirectR","PIndirectR")
  
  # Fit the mediation model
  med.fit <- (lm(M~X))
  out.fit <- (lm(Y~X+M))
  med.out <- mediate(med.fit,out.fit,treat = "X",mediator = "M",sims = nSimImai)
  
  # Get the direct and indirect effects
  mat_results["EDirectNR",1] <- summary(med.out)$z.avg
  mat_results["PDirectNR",1] <- summary(med.out)$z.avg.p
  mat_results["EIndirectNR",1] <- summary(med.out)$d.avg
  mat_results["PIndirectNR",1] <- summary(med.out)$d.avg.p
  
  # Reverse the order and run the mediation model again
  M2 <- Y
  Y2 <- M
  
  # Fit the mediation model
  med.fitR <- (lm(M2~X))
  out.fitR <- (lm(Y2~X+M2))
  med.outR <- mediate(med.fitR,out.fitR,treat = "X",mediator = "M2",sims = nSimImai)
  
  # Get the direct and indirect effects
  mat_results["EDirectR",1] <- summary(med.outR)$z.avg
  mat_results["PDirectR",1] <- summary(med.outR)$z.avg.p
  mat_results["EIndirectR",1] <- summary(med.outR)$d.avg
  mat_results["PIndirectR",1] <- summary(med.outR)$d.avg.p
  
  # Print out the matrix
  list(mat_results)
  
}
