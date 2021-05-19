normalizeK <- function(K,dX,dY){
 
  kVX = 1/sqrt(dX)
  kVY = 1/sqrt(dY)
  nK = K * (kVX %*% t(kVY))
  
  return(nK)
  
}
