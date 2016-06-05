#' @export

coxKernelnet <- function(x, y, t, nfolds, stdbeta, alpha)
{

# Correlation matrix
rbf <- as.matrix(getLaplacian(x, "RBF"))

# cross validation
foldid<-coxsplit(y, nfolds)

fit <- Coxnet(x, y, Omega = rbf, penalty="Net", alpha=alpha, foldid=foldid, isd=stdbeta)
beta <- fit$Beta

dt <- data.frame(cbind(y, x))
dt <- dt[order(dt[,1]),]
z <- dt[,3:ncol(dt)]

event <- subset(dt, dt$status==1 & dt$time <= t)

unique_time <- unique(event[,1])

R <- sapply(unique_time, function(t){
  which(dt[,1]>=t)})

h0 <- sum(sapply(1:length(unique_time), function(i){
  temp <- z[R[[i]],]
  1/sum((exp(as.matrix(temp) %*% as.matrix(beta))))
}))

s0 <- exp(-h0)

# Fetch the evaluation metrics
eval_metrics <- eval(x, y, t, beta, h0, s0, foldid, nfolds)

# Append the evaluation metrics to the Coxnet object and return
return(c(fit, s0=s0, unlist(eval_metrics)))

}