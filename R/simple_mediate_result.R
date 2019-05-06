setClass("SimpleMediateResult", 
         representation = representation(
           direct_p = "numeric",
           indirect_p="numeric"
         ),
         prototype = prototype(
           direct_p = 0,
           indirect_p = 0
         ),
         validity = function(object){
           if(object@direct_p<0 || object@direct_p>1){return ("direct p value must be between 0 and 1, inclusive")}
           if(object@indirect_p<0 || object@indirect_p>1){return ("indirect p value must be between 0 and 1, inclusive")}
         }
) -> SimpleMediateResult