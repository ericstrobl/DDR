CV_KRR_dotX_SQ3_fast <- function(dotX,Y,Hs,y_tr,cy){
  
  ## estimate conditional expectations
  n = nrow(Y); m = nrow(dotX)
  err = 99999999; #maximum error allowed
  for (s in seq_len(length(Hs))){
    for (l in seq_len(length(Hs[[1]]))){
      
      KRRs = t(t(Y) %*% Hs[[s]][[l]])
      res = LOO_KRR_fast(Y, Y-KRRs[1:n,],Hs[[s]][[l]])
      dens=rbind(res,KRRs[(n+1):m,]);
      
      dens = pmax(dens,0);
      dens = apply(dens, 1, function(x) x/trapz(cy,x)); dens = t(dens);
      term1t = apply(dens[1:m,],1,function(x) trapz(cy,x^2))
      # term1t = apply(dens[1:n,],1,function(x) cumsum(x^2))
      term1 = 0.5*mean(term1t)

      term2=c();
      for (e in 1:n){
        term2 =  c(term2, approx(x=cy,y=dens[e,],xout=y_tr[e],
                                 ties = "ordered")$y)
      }
      term2 = mean(term2);
      
      errn = term1 - term2;
      
      if (errn<err){
        err = errn
        pre = res
        pre_te = KRRs[(n+1):m,]
      }
    }
  }
  
  list(err=err,pre=pre,pre_te=pre_te);
  
}