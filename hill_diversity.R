#' Calculate Hill's true diversity metric
#'
#' @param abundances vector of species abundances
#' @param q this parameter defines the kind of mean used (Tuomisto 2012)
#'
#' @return single diversity value
#' @export
#'
#' @examples hill.diversity(c(1,20,34,32,355,5))
hill.diversity <- function(abundances,q = 1){
  abundances <- abundances[abundances != 0]
  abundances <- abundances/sum(abundances)
  R <- length(abundances)
  # hill diversity is not defined for q = 1,
  # but its limit exists and equals
  # the exponential of shannon entropy 
  # (Jost 2006,2007,Tuomisto 2012)
  if(q == 1){
    D <- exp(-sum(abundances*log(abundances)))
  }else{
    mean.p <- (sum(abundances*abundances^(q-1)))^(1/(q-1))
    D <- 1/mean.p
  }
  return(D)
}