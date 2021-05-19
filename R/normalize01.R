normalize01 <- function(X){
  X = as.matrix(X)
  
  for (d in 1:ncol(X)){
    X[,d] = (X[,d] - min(X[,d]))/(max(X[,d]) - min(X[,d]))
  }
  return(X)
}
