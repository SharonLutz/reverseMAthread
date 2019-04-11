suppressMessages(library(methods))

setClass("MediationProbValues",
         slots=representation(
           pval_direct="numeric",
           pval_indirect="numeric"
         ),
         prototype = prototype(
           pval_direct=0,
           pval_indirect=0
         ),
         validity = function(object){
           if(object@pval_direct<0 || object@pval_direct > 1){return("Error: pval_direct must be between 0.0 and 1.0")}
           if(object@pval_indirect<0 || object@pval_indirect > 1){return("Error: pval_indirect must be between 0.0 and 1.0")}
         })->MediationProbValues

setClass("MediateLinearModels",
            slots = representation(
              med.fit="lm",
              out.fit="lm"
            ),
            prototype=prototype(
              med.fit=NULL,
              out.fit=NULL
            )
            ) -> MediateLinearModels

setClass("MediateSimData",
         slots = representation(
           X = "numeric",
           M = "numeric",
           Y = "numeric",
           med.fit = "lm",
           out.fit = "lm"
         ),
         prototype = prototype(
           X = 0,
           M = 0,
           Y = 0,
           med.fit=NULL,
           out.fit=NULL
         )
         ) -> MediateSimData

setClass("MediateVariables", 
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
           if(object@n<0 | object@n==0 | floor(object@n)!=ceiling(object@n) ){return("Error: n must be an integer greater than or equal to 1")}
           if(object@nSim<0 | object@nSim==0 | floor(object@nSim)!=ceiling(object@nSim) ){return("Error: n must be an integer greater than or equal to 1")}
           
           #valid input checks
           
           if(object@pX<0 | object@pX>1){return("Error: pX must be greater than 0 and less than 1")}
           if(!object@varM>0){return("Error: varM must be greater than 0")}
           if(!object@varY>0){return("Error: varY must be greater than 0")}
           
           
           
           return(TRUE)
         }
         )->MediateVariables

setGeneric(name = "generateData",
           def = function(theObject, bM.ind){
             standardGeneric("generateData")
           })
setMethod(f="generateData",
          signature = "MediateVariables",
          definition = function(theObject, bM.ind=1){
            data = MediateSimData()
            
            data@X = rbinom(theObject@n,2,theObject@pX)
            
            data@M = rnorm(theObject@n,theObject@gamma0 + theObject@gammaX * data@X, sqrt(theObject@varM))
            
            data@Y = rnorm(theObject@n, theObject@beta0 + theObject@betaX * data@X + theObject@betaM[[bM.ind]] * data@M, sqrt(theObject@varY))
            
            return(data)
          })

setGeneric(name = "reverseData",
           def = function(theObject){
             standardGeneric("reverseData")
           })

setMethod(f="reverseData",
          signature = "MediateSimData",
          definition = function(theObject){
            tmpvar = theObject@M
            theObject@M = theObject@Y
            theObject@Y = tmpvar
            return(theObject)
          })

setGeneric(name = "generateDataMatrix",
           def = function(theObject){
             standardGeneric("generateDataMatrix")
           })

setMethod(f="generateDataMatrix",
          signature = "MediateVariables",
          definition = function(theObject){
            result = matrix(list(), nrow=length(theObject@betaM), ncol=theObject@nSim)
            
            for(i in 1:theObject@nSim){
              for(bM.ind in 1:length(theObject@betaM)){
                data_obj = generateData(theObject, bM.ind)
                result[[bM.ind,i]] = data_obj
              }
            }
            
            return(result)
          })

setGeneric(name = "reverseDataMatrix",
           def = function(theList){
             standardGeneric("reverseDataMatrix")
           })
setMethod(f="reverseDataMatrix",
          signature = "matrix",
          definition = function(theList){
            for(i in 1:length(theList)){
              theList[[i]] = reverseData(theList[[i]])
            }
            return(theList)
          })


setGeneric(name = "assembleLinearModels",
           def = function(theObject){
             standardGeneric("assembleLinearModels")
           })

setMethod(f="assembleLinearModels",
          signature = "MediateSimData",
          definition = function(theObject){
            df = data.frame(X=theObject@X, Y=theObject@Y, M=theObject@M)
            med.fit = lm("M~X", data=df)
            out.fit = lm("Y~X+M", data=df)
            return(MediateLinearModels(med.fit=med.fit, out.fit=out.fit))
          })
