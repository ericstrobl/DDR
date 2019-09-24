CV_KRR_dotX_SQ3_fast <- function(dotX,Y,Hs,cym,nn_cy){

  ## estimate conditional expectations
  n = length(Y); m = nrow(dotX)
  err = 99999999; #maximum error allowed

  hs = seq(0.1,0.5,length.out=20)

  est = foreach (h=hs,.export= c("LOO_KRR_fast", "trapz_fastM")) %dopar% {
  # est = for (h in hs){

    ty = exp( - ((Y - cym)^2) / (2*h^2) ) * (1/(sqrt(2*pi)*h))

    for (s in seq_len(length(Hs))){
      for (l in seq_len(length(Hs[[1]]))){

        KRRs = t(t(ty) %*% Hs[[s]][[l]])
        res = LOO_KRR_fast(ty, ty-KRRs[1:n,],Hs[[s]][[l]])
        dens = rbind(res,KRRs[(n+1):m,])

        dens = pmax(dens,0) + 1E-10;
        dens = dens/trapz_fastM(cym[1,2]-cym[1,1],dens)
        term1 = 0.5*mean(trapz_fastM(cym[1,2]-cym[1,1],dens^2))
        term2 = mean(dens[cbind(1:n,nn_cy)])
        errn = term1 - term2;

        if (errn<err){
          h_star = h
          err = errn
          pre = dens[1:n,]
          pre_te = dens[(n+1):m,]
        }
      }
    }
    list(err = err, pre = pre, pre_te = pre_te, h_star=h_star)
  }

  err = 9999999;
  for (e in seq_len(length(est))){
    if (est[[e]]$err < err){
      densf=rbind(est[[e]]$pre, est[[e]]$pre_te)
      h_star = est[[e]]$h_star
      err = est[[e]]$err
    }
  }

  return(list(densf=densf,err=err,h_star=h_star))

}
