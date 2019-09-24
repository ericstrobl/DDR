#' DDR - estimates conditional densities using Gaussian kernel
#' @param X_tr training predictors matrix (samples by features)
#' @param y_tr training response vector
#' @param X_te test predictors matrix (samples by features)
#' @return A list containing the response values \code{y} and
#' their conditional density estimates \code{dens} for all samples in the test set
#' @export
#' @examples
#' DDR(rnorm(100),rnorm(100),rnorm(100));
#'

DDR <- function(X_tr,y_tr,X_te){
  ncy = 500;

  X_tr = t(t(X_tr))
  X_te = t(t(X_te))

  n=nrow(X_tr)
  ty = matrix(NaN,n,ncy)
  dens = rep(NaN,ncy)

  #normalize predictors
  X = rbind(X_tr,X_te);
  X = normalize(X)
  m=nrow(X);

  #normalize response
  my = mean(y_tr);
  sy = sd(y_tr);
  y_tr = (y_tr - my)/sy
  cy = seq(min(y_tr),max(y_tr),length.out=ncy);

  nn_cy = get.knnx(cy,y_tr,k=1)$nn.index #for fast computation of term 2
  cym = matrix(cy); cym = repmat(t(cy),n,1);

  #compute RBF kernels
  dotX=(dist(cbind(X,X),diag = TRUE, upper = TRUE))^2; dotX = as.matrix(dotX);
  med_dist = median(c(t(dist(X))))
  sigmas = med_dist*seq(0.5,2,0.3);
  lambdas = c(1E-1,1E-2,1E-3,1E-4,1E-5,1E-6);
  Hs = lapply(1:length(sigmas),function(x) vector("list", length = length(lambdas)));
  for (s in seq_len(length(sigmas))){
    K = RBF_kernel(dotX,sigmas[s])
    for (l in seq_len(length(lambdas))){
      K_tr = K[1:n,1:n]
      Hs[[s]][[l]] = chol2inv(chol( K_tr + m*diag(nrow(K_tr))*lambdas[l])) %*%
        K[1:n,1:m];
    }
  }

  # perform DDR in parallel (Step 1)
  out=CV_KRR_dotX_SQ3_fast(dotX,y_tr,Hs,cym,nn_cy)
  h_star = out$h_star

  # sharpen estimate (Step 3)
  densf = sharpen_density(out$densf,cym,nn_cy)
  dens = densf[(n+1):m,]

  cy = sy*cy + my;
  dens = dens / sy

  return(list(y = cy, dens = dens))

}
