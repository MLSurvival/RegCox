eval <- function(x, y, t, beta, h0, s0, foldid, nfolds)
{
  
  x_size <- dim(x)[1]
  y_size <- dim(x)[2]
  events <- sum(y[,2])
  
  # Find the prediction on the training data
  ht <- -log(s0^(exp(x%*%as.matrix(beta))))
  st <- s0^(exp(x%*%as.matrix(beta)))
  p <- as.list(ht)

  # AUC calculation
  auc <- sapply(c(1:nfolds), function(n){
    x1 <- ht[which(foldid==n),]
    t1 <- y[which(foldid==n),1]
    e1 <- y[which(foldid==n),2]
    output <- concordance.index(x=x1, surv.time=t1, surv.event=e1, method="noether")
    return(output$c.index)
  })
  
  # R square calculation
  R2 <- 1-exp((2/events)*(coxloglik(x_size, y, x, rep(0, y_size)) - coxloglik(x_size, y, x, beta)))

  # Martingale Residual calculation
  MR <- sapply(c(1:nfolds), function(n) { 
    sum((y[which(foldid==n),2] - (exp(x[which(foldid==n), ] %*% beta)*h0))^2
    )/sum(foldid==n) } )
  
  brier_score = sum((y[,2]-st)^2)/x_size
  
  return(c(auc=mean(auc), auc_sd=sd(auc), R2=R2, MR=mean(MR), MR_sd = sd(MR), brier_score = brier_score[[1]]))
}