normalize <- function(mat){
  
  if (is.null(nrow(mat))){mat = matrix(mat);}
  
  mat = apply(mat, 2, function(x) if (sd(x)>0){(x - mean(x)) / sd(x)} else{x-mean(x);})
  
  
}