
# R implementation of C++ utility function used in umod.cpp

ufun_labouR_disc <- function(e,s,l,params){
  if (is.null(dim(e))) e <- matrix(e,length(e),1)
  s <- matrix(s,dim(e)[1],1)
  # return objects
  util <- matrix(0,dim(e)[1],dim(e)[2])
  phivals <- c(0,params$phival,1)
  phivec <- rep(0,length(s))
  for (i in 1:length(s)){
    phivec[i] <- phivals[ s[i] + 1 ]
  }
  phimat <- matrix(phivec,length(s),ncol(e))
  hfac <- exp( params$theta * phimat )
  for (i in 1:nrow(e)){
    for (j in 1:ncol(e)){
      dist <- e[i,j] - params$cutoff
      if (dist>=0){
        # positive consumption?
        util[i,j] <- hfac[i,j] * (1/(1-params$gamma)) * (e[i,j] * exp(params$alpha*l))^(1-params$gamma)  
      } else {
        z <- hfac[i,j] *exp(l * params$alpha)^(1-params$gamma) 
        grad      <-  params$cutoff^(-params$gamma) * z
        hess      <- (-params$gamma/params$cutoff) * grad
        util[i,j] <- (1/(1-params$gamma)) * params$cutoff^(1-params$gamma) * z + grad*dist + 0.5 * hess * dist^2
      }
    }
  }
  util <- util + params$mu * phimat
  return(util)
}
