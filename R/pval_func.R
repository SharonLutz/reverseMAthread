pval <- function(x, xhat){
  ## Compute p-values
  if (xhat == 0) out <- 1
  else {
    out <- 2 * min(sum(x > 0), sum(x < 0)) / length(x)
  }
  return(min(out, 1))
}
