LOO_KRR_fast <- function(tr_y, tr_e, H){
  err = tr_e/(1-diag(H))
  pre = -(err - tr_y)
  
  pre
}