#' @export

cocktail.cmd <- function(x, y, t, nfolds, stdbeta, tol, lambda, itermax, alpha, weightj){

if(!stdbeta)
{
	print("ERROR: stdbeta must be TRUE for Cocktail regularization")
	return()
}

beta.old <- rep(0, ncol(x))

D <- cocktail.doj(nrow(x), x, y)

for(i in 1:itermax){

G <- cocktail.gradient(nrow(x), y, x, beta.old)
beta.new <- sapply(c(1:ncol(x)), function(j){
  soft.threshold(D[j]*beta.old[j]-G[j], lambda*alpha*weightj)/
            ( D[j] + lambda*(1-alpha) )
})

beta.new <- (beta.new-min(beta.new))/(max(beta.new)-min(beta.new))

diff <- beta.new - beta.old
val <- norm(as.matrix(diff),"o")

if (val < tol)
  break


beta.old <- beta.new

}

beta <- beta.new

	foldid <- coxsplit(y, nfolds)
  
	dt <- data.frame(cbind(y, x))
	dt <- dt[order(dt[,1]),]
	xtemp <- dt[,3:ncol(dt)]
	
	event <- subset(dt, dt$status==1 & dt$time <= t)
	
	unique_time <- unique(event[,1])
	
	R <- sapply(unique_time, function(t){
	  which(dt[,1]>=t)})
	
	h0 <- sum(sapply(1:length(unique_time), function(i){
	  temp <- xtemp[R[[i]],]
	  1/sum((exp(as.matrix(temp) %*% as.matrix(beta))))
	}))
	
	s0 <- exp(-h0)
	
	# Fetch the evaluation metrics
	eval_metrics <- eval(x, y, t, beta, h0, s0, foldid, nfolds)
	
	# Append the evaluation metrics to the Coxnet object and return
	return(c(list(Beta = beta), s0=s0, unlist(eval_metrics)))

}


cocktail.gradient = function(n,y,z,beta){
  
  # combine x and y and pick first occurence of x for unique time in y
  dt <- data.frame(cbind(y,z))
  dt <- dt[order(dt[,1]),]
  z <- dt[,3:ncol(dt)]
  eventx <- subset(dt, dt$status==1)
  
  unique_time <- unique(eventx[,1])
  
  R <- sapply(unique(eventx[,1]), function(t){
    which(dt[,1]>=t)
  })
  
  grad <- sapply(c(1:ncol(z)), function(j){
    
    sum(sapply(c(1:length(unique_time)), function(i){
      
      p1 <- subset(eventx, eventx$time==unique_time[i])[1,j]
      
      p2 <- (t(as.matrix(z[R[[i]],j]))) %*% (exp(as.matrix(z[R[[i]],]) %*% as.matrix(beta)))
      
      p3 <- sum(exp(as.matrix(z[R[[i]],]) %*% as.matrix(beta)))
      
      return((p2/p3)-p1)
      
    }))/n
    
  })
  
  return(grad)

}


cocktail.doj <- function(n, x, y)
{
  
  # combine x and y and pick first occurence of x for unique time in y
  dt <- data.frame(cbind(y,x))
  dt <- dt[order(dt[,1]),]
  z <- dt[,3:ncol(dt)]
  eventx <- subset(dt, dt$status==1)
  
  unique_time <- unique(eventx[,1])
  
  R <- sapply(unique(eventx[,1]), function(t){
    which(dt[,1]>=t)
  })
  
  D <- sapply(c(1:ncol(z)), function(j){
    
    sum(sapply(c(1:length(unique_time)), function(i){
      
      p1 <- max(x[R[[i]],j])
      
      p2 <- min(x[R[[i]],j])
      
      return((p1-p2)^2)
      
    }))/(4*n)
  })

  return(D)
}


soft.threshold<-function(vec,lam)
  {
  if ( length(lam)>1 & length(lam)!=length(vec) ) {
    cat('\n ERROR: THE SIZE OF THE SECOND ARGUMENT SHOULD BE 1 OR THE SAME AS THE SIZE OF THE FIRST ARGUMENT.\n')
    return ( 0 )
  }
  idx.1<-which(vec < -lam)
  idx.2<-which(vec > lam)
  res<-rep(0,length(vec))
  if (length(lam)==1){
    if ( length(idx.1)>0 ) res[idx.1]<- vec[idx.1]+lam
    if ( length(idx.2)>0 ) res[idx.2]<- vec[idx.2]-lam
  } else if (length(lam)>1) {
    if ( length(idx.1)>0 ) res[idx.1]<- vec[idx.1]+lam[idx.1]
    if ( length(idx.2)>0 ) res[idx.2]<- vec[idx.2]-lam[idx.2]
  }
  return( res )
}

