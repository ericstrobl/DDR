sharpen_density <- function(dens,cym,nn_cy){

  err=9999999;
  n=length(nn_cy)
  delta_grid = seq(0, 0.5, length.out = 50)
  for (d in 1:50){
    denst = dens;
    denst = pmax(denst - delta_grid[d],0)+1E-10;
    denst = denst/trapz_fastM(cym[1,2]-cym[1,1],denst)

    term1 = 0.5*mean(trapz_fastM(cym[1,2]-cym[1,1],denst^2))

    if (!(is.nan(term1))){
      term2 = mean(denst[cbind(1:n,nn_cy)])
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

  return(densf)

}
