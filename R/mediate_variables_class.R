library(methods)

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
#n=6000
#pX=0.2
#gamma0=0
#gammaX=0.1
#varM=1
#beta0=0
#betaX=1
#betaM=c(0,0.1,0.2)
#varY=1
#nSim=1000
#nSimImai=10000
#SEED=1
#alpha_level=0.05
#bM.ind <- 1

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
           nSim = "integer",
           nSimImai = "integer",
           SEED = "integer",
           alpha_level = "numeric",
           bM.ind = "numeric"
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
           nSim = as.integer(1000),
           nSimImai = as.integer(10000),
           SEED = as.integer(1),
           alpha_level=0.05,
           bM.ind = 1
         ),
         validity = function(object){
           if(length(object@n) != 1){return ("slot n must be a single numeric value")}
           if(length(object@pX) != 1){return ("slot pX must be a single numeric value")}
           if(length(object@gamma0) != 1){return ("slot gamma0 must be a single numeric value")}
           if(length(object@gammaX) != 1){return ("slot gammaX must be a single numeric value")}
           if(length(object@varM) != 1){return ("slot varM must be a single numeric value")}
           if(length(object@beta0) != 1){return ("slot beta0 must be a single numeric value")}
           if(length(object@betaX) != 1){return ("slot betaX must be a single numeric value")}
           if(length(object@betaM) != 3){return ("slot betaM must be a vector of 3 numeric values")}
           if(length(object@varY) != 1){return ("slot varY must be a single numeric value")}
           if(length(object@nSim) != 1){return ("slot nSim must be a single integer value")}
           if(length(object@nSimImai) != 1){return ("slot nSimImai must be a single integer value")}
           if(length(object@SEED) != 1){return ("slot SEED must be a single integer value")}
           if(length(object@alpha_level) != 1){return ("slot alpha_level must be a single numeric value")}
           if(length(object@bM.ind) != 1){return ("slot bM.ind must be a single numeric value")}
           return(TRUE)
         }
         )->MediateVariables

setGeneric(name = "generateData",
           def = function(theObject){
             standardGeneric("generateData")
           })
setMethod(f="generateData",
          signature = "MediateVariables",
          definition = function(theObject){
            data = MediateSimData()
            set.seed(theObject@SEED)
            
            data@X = rbinom(
              theObject@n,
              2,
              theObject@pX
              )
            
            data@M = rnorm(
              theObject@n,
              theObject@gamma0 + theObject@gammaX * data@X, 
              sqrt(theObject@varM)
              )
            
            data@Y = rnorm(
              theObject@n, 
              theObject@beta0 + theObject@betaX * data@X + theObject@betaM[theObject@bM.ind] * data@M, 
              sqrt(theObject@varY)
              )
            return(data)
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
