#' Active learning using regularized cox models
#' @param filename Input matrix. Each row is an observation vector.
#' @param regtype Type of regularization to apply {"Lasso", "Enet", "Lapnet", "Kernelnet"}
#' @param metric Type of evaluation metric to fetch {"AUC", "MSE"}
#' @param nfolds Number of folds for cross validation
#' @param standardize Logical flag for standardizing x using z-score normalization. default = TRUE
#' @param alpha ratio between L1 and Laplacian for "Lapnet", or between L1 and L2 for "Enet". Not needed for "Lasso". Value should be between 0 and 1
#' @param lambda a user supplied decreasing sequence. can range from 0.01 to 1
#' @export

arc <- function(filename, regtype, metric, nfolds, standardize=TRUE, alpha, lambda)
{

data <- read.csv(filename, head=TRUE, sep=",")

if(standardize){
x <- as.matrix(data[,3:ncol(data)])
y <- as.matrix(data[,1:2])
xnorm <- apply(x, 2, scale)
data <- data.frame(cbind(y, xnorm))
}

if(metric == "AUC"){
  if(regtype == "Lasso")
    return(arc.lasso.auc(data, nfolds, lambda))
  if(regtype == "Enet")
    return(arc.enet.auc(data, nfolds, alpha, lambda))
  if(regtype == "Lapnet")
    return(arc.lapnet.auc(data, nfolds, alpha, lambda))
  if(regtype == "Kernelnet")
    return(arc.kernelnet.auc(data, nfolds, alpha, lambda))
  else
    print("Unknown regularization type")}
if(metric == "MSE"){
  if(regtype == "Lasso")
    return(arc.lasso.mse(data, nfolds, lambda))
  if(regtype == "Enet")
    return(arc.enet.mse(data, nfolds, alpha, lambda))
  if(regtype == "Lapnet")
    return(arc.lapnet.mse(data, nfolds, alpha, lambda))
  if(regtype == "Kernelnet")
    return(arc.kernelnet.mse(data, nfolds, alpha, lambda))
  else
    print("Unknown regularization type")}
else
    print("Unknown metric")

}
