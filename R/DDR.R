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
      Hs[[s]][[l]] = chol2inv(chol( K_tr + diag(nrow(K_tr))*lambdas[l])) %*%
        K[1:n,1:m];
    }
  }

  # perform DDR in parallel (Step 1)
  hs = seq(0.1,0.5,length.out=20)

  dens = matrix(NaN,m,ncy)
  densf = dens[(n+1):m,];
  est = foreach (h=hs,.export= c("CV_KRR_dotX_SQ3_fast", "LOO_KRR_fast"),
                 .packages=c("pracma")) %dopar% {
                   # est = for(h in hs){
                   #   print(h)

                   ty = exp( - ((y_tr - cym)^2) / (2*h^2) ) * (1/(sqrt(2*pi)*h))

                   out=CV_KRR_dotX_SQ3_fast(dotX,ty,Hs,y_tr,cy)

                   dens[1:n,]=out$pre;
                   dens[(n+1):m,]=out$pre_te

                   list(densf=dens,loss_pre=out$err)
                 }

  loss_pre = 9999999;
  for (e in seq_len(length(est))){
    if (est[[e]]$loss_pre < loss_pre){
      densf=est[[e]]$dens;
      loss_pre = est[[e]]$loss_pre
    }
  }

  # enforce non-negativity and AUC of 1 (Step 2)
  densf = pmax(densf,0)
  dens=densf;

  dens = apply(dens, 1, function (x) x/trapz(cy,x));
  dens=t(dens);

  # sharpen estimate (Step 3)
  err=9999999;
  delta_grid = seq(0, 0.5, length.out = 50)
  for (d in 1:50){
    denst = dens;
    denst = pmax(denst - delta_grid[d],0)
    denst = apply(denst,1,function(x) x/trapz(cy,x));
    denst = t(denst)

    term1t = apply(denst[1:m,],1,function(x) trapz(cy,x^2))
    # term1t = apply(dens[1:n,],1,function(x) cumsum(x^2))
    term1 = 0.5*mean(term1t)

    term2=c();
    if (!(is.nan(term1))){
      for (e in 1:n){
        term2 =  c(term2, approx(x=cy,y=denst[e,],xout=y_tr[e],
                                 ties = "ordered")$y)
      }
      term2 = mean(term2);
    } else{
      term2 = NaN
    }


    errn = term1 - term2;
    if (is.nan(errn)){
      break
    }
    if (errn < err){
      densf = denst
      err = errn;
    }
  }
  dens = densf[(n+1):m,];



  cy = sy*cy + my;

  dens = dens / sy
  # dens = pmin(pmax(0,dens),1/(sqrt(2*pi)*h));
  # plot(cy,dnorm(cy,1,sqrt(0.2)),col='blue',ylim=c(0,max(dens)+0.001))
  # lines(cy,dens,xlim=c(min(y_tr),max(y_tr)))

  return(list(y = cy, dens = dens))

}
