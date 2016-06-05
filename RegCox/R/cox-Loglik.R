coxloglik = function(n,y,z,beta){
  
  delta <- y[,2]
  
  L <- delta%*%(z%*%as.matrix(beta)) -
    
    delta%*%log(sapply(y[,1], function(p){
      # find instances of R
      # foreach j in R
      R <- z[which(y[,1] >= p),]
      sum(exp(R %*% beta))
    }))
  
  return(-L/n)
}