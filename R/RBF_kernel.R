RBF_kernel <- function(dotx,sigma){
  
  kx = exp(-dotx/(2*sigma^2));
  
  return(kx)
  
}