
#' @title MediateDataGenerationParameters
#' @description Convenience class for use with other methods and classes to generate mediation data for use with parallel processing.
#' @author Michael Gooch
#' @name MediateDataGenerationParameters
#' @slot n is the sample size.
#' @slot pX is the minor allele frequency
#' @slot gamma0 is the intercept for M
#' @slot gammaX is the association of X with M
#' @slot varM is the variance of M
#' @slot beta0 is the intercept for Y
#' @slot betaX is the direct effect of X on Y
#' @slot betaM is a vector of different associations of M with Y
#' @slot varY is the variance of Y
#' @slot nSim is the number of simulations to run
methods::setClass("MediateDataGenerationParameters", 
         slots = representation(
           n = "numeric",
           pX = "numeric",
           gamma0 = "numeric",
           gammaX = "numeric",
           varM = "numeric",
           beta0 = "numeric",
           betaX = "numeric",
           betaM = "numeric",
           varY = "numeric",
           nSim = "numeric"
         ),
         prototype = prototype(
           n = 6000,
           pX = 0.2,
           gamma0 = 0,
           gammaX = 0.1,
           varM = 1,
           beta0 = 0,
           betaX = 1,
           betaM = c(0,0.1,0.2),
           varY = 1,
           nSim = 1000
         ),
         validity = function(object){
           # Error checks
           
           #number of values checks
           if(length(object@n) != 1){return ("Error:  n must be a single numeric value")}
           if(length(object@pX) != 1){return ("Error: pX must be a single numeric value")}
           if(length(object@gamma0) != 1){return ("Error: gamma0 must be a single numeric value")}
           if(length(object@gammaX) != 1){return ("Error: gammaX must be a single numeric value")}
           if(length(object@varM) != 1){return ("Error: varM must be a single numeric value")}
           if(length(object@beta0) != 1){return ("Error: beta0 must be a single numeric value")}
           if(length(object@betaX) != 1){return ("Error: betaX must be a single numeric value")}
           if(length(object@varY) != 1){return ("Error: varY must be a single numeric value")}
           if(length(object@nSim) != 1){return ("Error: nSim must be a single integer value")}
           
           
           if(length(unique(object@betaM)) != length((object@betaM))){return("Error: values in betaM must be unique")}
           if(length(unique(object@betaM))<2){return("Error: betaM must be a vector with at least two values")}
           
           #INT style values
           if(object@n<=0 || floor(object@n)!=ceiling(object@n) ){return("Error: n must be an integer greater than or equal to 1")}
           if(object@nSim<=0 || floor(object@nSim)!=ceiling(object@nSim) ){return("Error: n must be an integer greater than or equal to 1")}
           
           #valid input checks
           
           if(object@pX<0 | object@pX>1){return("Error: pX must be greater than 0 and less than 1")}
           if(object@varM<=0){return("Error: varM must be greater than 0")}
           if(object@varY<=0){return("Error: varY must be greater than 0")}
           
           return(TRUE)
         })->MediateDataGenerationParameters

#' @title MediateModelVariables
#' @description convenience container for vectors of simulated data to be used with other methods and classes to generate linear models for use with mediation analysis in parallel processes
#' @author Michael Gooch
#' @name MediateModelVariables
#' @slot X is the exposure
#' @slot M is the mediator
#' @slot Y is the outcome
#' @slot SEED is the seed to use when performing mediation on the models generated from this data
methods::setClass("MediateModelVariables",
         slots = representation(
           X = "numeric",
           M = "numeric",
           Y = "numeric",
           SEED = "numeric"
         ),
         prototype = prototype(
           X = 0,
           M = 0,
           Y = 0,
           SEED = 1
         )
) -> MediateModelVariables

#' @title MediateLinearModels
#' @description convenience container for 4 linear models, 2 pairs, intended for forward and reverse mediation analysis
#' @author Michael Gooch
#' @name MediateLinearModels
#' @slot med.fit linear model for the mediator
#' @slot out.fit linear model for the outcome
#' @slot med.fit.r linear model for the mediator with M and Y swapped
#' @slot out.fit.r linear model for the outcome with M and Y swapped
methods::setClass("MediateLinearModels",
         slots = representation(
           med.fit="lm",
           out.fit="lm",
           med.fit.r="lm",
           out.fit.r="lm"
         ),
         prototype=prototype(
           med.fit=NULL,
           out.fit=NULL,
           med.fit.r=NULL,
           out.fit.r=NULL
         )
) -> MediateLinearModels

#' @title MediationProbValues
#' @description Convenience container for the result of a forward and reverse pair of mediation analysis results
#' @author Michael Gooch
#' @name MediationProbValues
#' @slot pval_direct direct probability from forward mediation analysis
#' @slot pval_indirect indirect probability from forward mediation analysis
#' @slot pval_direct_r direct probability from reverse mediation analysis
#' @slot pval_indirect_r indirect probability from reverse mediation analysis
methods::setClass("MediationProbValues",
         slots=representation(
           pval_direct="numeric",
           pval_indirect="numeric",
           pval_direct_r="numeric",
           pval_indirect_r="numeric"
         ),
         prototype = prototype(
           pval_direct=0,
           pval_indirect=0,
           pval_direct_r=0,
           pval_indirect_r=0
         ),
         validity = function(object){
           if(object@pval_direct<0 || object@pval_direct > 1){return("Error: pval_direct must be between 0.0 and 1.0")}
           if(object@pval_indirect<0 || object@pval_indirect > 1){return("Error: pval_indirect must be between 0.0 and 1.0")}
           if(object@pval_direct_r<0 || object@pval_direct_r > 1){return("Error: pval_direct_r must be between 0.0 and 1.0")}
           if(object@pval_indirect_r<0 || object@pval_indirect_r > 1){return("Error: pval_indirect_r must be between 0.0 and 1.0")}
         })->MediationProbValues

#' @name generateData
#' @title generateData
#' @author Michael Gooch
#' @description use data generation parameters with the random number generator to generate simulated data
#' @param theObject the object to be used in generation of the data, should be an instance of \code{MediateDataGenerationParameters}
#' @param bM.ind index for the betaM vector stored in theObject
#' @param SEED SEED to use with the resulting data
#' @return an instance of \code{MediateModelVariables} with random values determined by current rng state, MediateDataGenerationParameters, and bM.ind,
methods::setGeneric(name = "generateData",
           def = function(theObject, bM.ind=1, SEED=1){
             standardGeneric("generateData")
           })

methods::setMethod(f="generateData",
          signature = "MediateDataGenerationParameters",
          definition = function(theObject, bM.ind=1, SEED=1){
            
            X = rbinom(theObject@n,2,theObject@pX)
            
            M = rnorm(theObject@n,theObject@gamma0 + theObject@gammaX * X, sqrt(theObject@varM))
            
            Y = rnorm(theObject@n, theObject@beta0 + theObject@betaX * X + theObject@betaM[[bM.ind]] * M, sqrt(theObject@varY))
            
            return(MediateModelVariables(X=X, M=M, Y=Y, SEED=SEED))
          })

#' @name generateDataMatrix
#' @title generateDataMatrix
#' @author Michael Gooch
#' @description generate all the data vectors and rng seeds needed to do mediation with parallel processing 
#' @param theObject an instance of MediateDataGenerationParameters
#' @param initial_SEED the SEED that will be used with set.seed before generating the data, the generated MediateModelVariables will be provided incremented SEED values starting with initial_SEED + 1
#' @return matrix with dimensions of \code{rows = length(theObject@@betaM)} and \code{columns = theObject@@nSim} containing the generated data
methods::setGeneric(name = "generateDataMatrix",
           def = function(theObject, initial_SEED=1){
             standardGeneric("generateDataMatrix")
           })

methods::setMethod(f="generateDataMatrix",
          signature = "MediateDataGenerationParameters",
          definition = function(theObject, initial_SEED=1){
            result = matrix(list(), nrow=length(theObject@betaM), ncol=theObject@nSim)
            set.seed(initial_SEED)
            next_seed = initial_SEED + 1
            for(i in 1:theObject@nSim){
              for(bM.ind in 1:length(theObject@betaM)){
                result[[bM.ind,i]] = generateData(theObject, bM.ind, next_seed)
                next_seed = next_seed + 1
              }
            }
            
            return(result)
          })

#' @name assembleLinearModels
#' @title assembleLinearModels
#' @author Michael Gooch
#' @description create forward and reverse(Y and M are swapped) linear models for use with mediation analysis
#' @param theObject an instance of MediateModelVariables
#' @return an instance of MediateLinearModels with forward and reverse linear models created using the variable vectors in theObject
methods::setGeneric(name = "assembleLinearModels",
           def = function(theObject){
             standardGeneric("assembleLinearModels")
           })

methods::setMethod(f="assembleLinearModels",
          signature = "MediateModelVariables",
          definition = function(theObject){
            df = data.frame(X=theObject@X, Y1=theObject@Y, M1=theObject@M, M2=theObject@Y, Y2=theObject@M)
            med.fit = lm("M1~X", data=df)
            out.fit = lm("Y1~X+M1", data=df)
            med.fit.r = lm("M2~X", data=df)
            out.fit.r = lm("Y2~X+M2", data=df)
            return(MediateLinearModels(med.fit=med.fit, out.fit=out.fit, med.fit.r=med.fit.r, out.fit.r=out.fit.r))
          })
