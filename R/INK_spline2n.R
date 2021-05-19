INK_spline2n <- function(X,Y){
  
  X = as.matrix(X)
  Y = as.matrix(Y)
  
  nX = nrow(X)
  nY = nrow(Y)
  cX = ncol(X)
  
  dX = matrix(0,nX,cX)
  for (x in 1:cX){
    dX[,x] = (1/3)*X[,x,drop=FALSE]^3 + X[,x,drop=FALSE]^2 + 1
  }
  
  dY = matrix(0,nY,cX)
  for (y in 1:cX){
    dY[,y] = (1/3)*Y[,y,drop=FALSE]^3 + Y[,y,drop=FALSE]^2 + 1
  }
  
  K = matrix(1,nX,nY)
  for (x in 1:cX){
    mm = outer(X[,x],Y[,x],FUN="pmin")
    Kt = (1/3)*mm^3 + 0.5*abs(outer(X[,x],Y[,x],FUN="-"))*mm^2 + 
      X[,x,drop=FALSE] %*% t(Y[,x,drop=FALSE]) + 1

    K = K * normalizeK(Kt,dX[,x,drop=FALSE],dY[,x,drop=FALSE])
  }
  
  return(K)
  
}
