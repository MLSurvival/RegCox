getLaplacian <- function(x, metric="LAPLACIAN")
{
  
  # RBF kernel calculation
  S <- sapply(1:dim(x)[2], function(m){
    sapply(1:dim(x)[2], function(n){
      exp(-1*(sqrt(sum((x[,m]-x[,n])^2))/(2*.09))) }) })
  diag(S) <- 0
  
  if(metric=="RBF")
  {
    return(S)
  }
  
  # Laplacian Calculation
  D <- diag(colSums(S))
  delta <- .5
  I <- diag(dim(x)[2])
  D[is.infinite(D^-.5)] <- 0
  Laplacian <- (1+delta)*I-(D%*%S%*%D)
  return(Laplacian)
  
}