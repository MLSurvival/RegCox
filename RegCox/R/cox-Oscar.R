#' @export

regcoxfista=function(delta, z, time, t, y_target, l1, l2, nfolds, stdbeta){
        beta=numeric(ncol(z))
        t1=1
        niter=6 # number of iterations
        y<-beta
	L<-100	# setting L value
	lambda1 <- l1 # L1 regularization parameter
	lambda2 <- l2 # pairwise L infty regularization parameter
	for(i in 1:niter){
        g <- dloglik(nrow(z),delta,z,y,time)
        v=y-(g/L)
        prox_flat <-  oscarapo(v,lambda1/L,lambda2/L) # oscar approximate proximal operator
        betanew <- prox_flat
        tnew<- (1 + sqrt(1 + 4*t1^2))/2
        y = as.vector(betanew + ((t1-1)/(tnew))*(betanew-beta))
        beta=betanew
        t1 = tnew;
        }
  
	if(stdbeta)
		beta <- (beta-min(beta))/(max(beta)-min(beta))
		
	foldid <- coxsplit(y_target, nfolds)
  
	dt <- data.frame(cbind(y_target, z))
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
	eval_metrics <- eval(z, y_target, t, beta, h0, s0, foldid, nfolds)
	
	# Append the evaluation metrics to the Coxnet object and return
	return(c(list(Beta = beta), s0=s0, unlist(eval_metrics)))
	
}

#gradient function
dloglik = function(n,delta,z,beta,time){
term1 <- t(delta)%*%z
  term2_p <- sapply(time, function(p){
      R <- z[which(time >= p),]
      if(!is.null(nrow(R))){
      t(R)%*%exp(R %*% beta)/sum(exp(R %*% beta))
      }
      else{
      as.matrix(numeric(length(beta)))
      }
    })
  term2 <- t(delta)%*%t(term2_p)
  L <- -1*(term1-term2)
return(L/n) 
}  

#function to create a permutation
permutation <- function(x){
absolute_x <- abs(x)
x1 <- sort(absolute_x,index.return=TRUE,decreasing=TRUE)
ordr_indx <- x1$ix
p <- matrix(0,nrow=length(ordr_indx),ncol=length(ordr_indx))
for(i in 1:length(ordr_indx)){
p[ordr_indx[i],i] <- 1
}
p_x <- list(x1$x,p)
names(p_x) <- c('x_tilda','permutation')
p_x[[1]] <- x%*%p
return(p_x)
}

#Approximate proximal operator for OSCAR
oscarapo <- function(v,lambda1,lambda2){
n <- length(v)
w <- lambda1+lambda2*(n-as.matrix(c(1:n)))
p_x <- permutation(v)
term1_u <- sin(p_x$x_tilda)
vec_1 <- as.vector(t(p_x$x_tilda)-w)
j <- length(vec_1)
term2_u <- numeric()
for(i in 1:j){
term2_u <- c(term2_u,max(vec_1[i],0))
}
u <- term1_u * t(term2_u)
X_star <- t(p_x$permutation)%*%t(u)
X_star
}

normalize = function(x){
  y = (x-mean(x))/sqrt(sum((x-mean(x))^2))
  return(y)
}

