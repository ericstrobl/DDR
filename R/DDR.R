#' DDR - estimates conditional densities using Gaussian kernel
#' @param X_tr training predictors matrix (samples by features)
#' @param y_tr training response vector
#' @param X_te test predictors matrix (samples by features)
#' @param lb lower bound of response variable
#' @param ub upper bound of response variable
#' @return A list containing the response values \code{y} and
#' their conditional density estimates \code{dens} for all samples in the test set
#' @export
#' @examples
#' DDR(rnorm(100),rnorm(100),rnorm(100));
#'

DDR <- function(X_tr,y_tr,X_te,lb=-Inf,ub=Inf){
  ncy = 500;
  
  X_tr = t(t(X_tr))
  X_te = t(t(X_te))
  
  n=nrow(X_tr)
  ty = matrix(NaN,n,ncy)
  dens = rep(NaN,ncy)
  
  #normalize predictors
  X = rbind(X_tr,X_te);
  X = normalize01(X)
  m = nrow(X);
  
  #normalize response
  my = mean(y_tr);
  sy = sd(y_tr);
  y_tr = (y_tr - my)/sy
  lb = (lb - my)/sy
  ub = (ub - my)/sy
  cy = seq(min(y_tr),max(y_tr),length.out=ncy);
  
  nn_cy = get.knnx(cy,y_tr,k=1)$nn.index #for fast computation of term 2
  cym = matrix(cy); cym = repmat(t(cy),n,1);
  
  #compute RBF kernels
  lambdas = c(1E-1,1E-2,1E-3,1E-4,1E-5,1E-6);
  # lambdas = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  Hs = lapply(1,function(x) vector("list", length = length(lambdas)));
  K = INK_spline2n(X[1:n,,drop=FALSE],X)
  for (l in seq_len(length(lambdas))){
    K_tr = K[1:n,1:n]
    
    Hs[[1]][[l]] = spdinv( K_tr + m*diag(nrow(K_tr))*lambdas[l]) %*%
    K[1:n,1:m];
    
    # Hs[[1]][[l]] = spdinv( (1-lambdas[l])*K_tr + diag(nrow(K_tr))*lambdas[l]) %*%
    # K[1:n,1:m];
  }
  
  # perform DDR in parallel (Step 1)
  out=CV_KRR_dotX_SQ3_fast_bound(X,y_tr,Hs,cym,nn_cy,lb,ub)
  h_star = out$h_star
  
  # sharpen estimate (Step 3)
  densf = sharpen_density(out$densf,cym,nn_cy)
  dens = densf[(n+1):m,]
  
  cy = sy*cy + my;
  dens = dens / sy
  
  return(list(y = cy, dens = dens))
  
}
