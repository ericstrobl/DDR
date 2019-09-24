trapz_fastM <- function(delta_x, Y){
  ly = ncol(Y);
  return(0.5*delta_x*(Y[,1] + Y[,ly]+ 2 * rowSums(Y[,2:(ly-1)])))
}
