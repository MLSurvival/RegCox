install_github("dmkdcs/RegCox/RegCox")
library(RegCox)

# Set the path to the file below
filename = "synthetic1-survival.csv"
#filename = "breast.csv"

train_time <- train_data$time
train_event <-train_data$event
x <- as.matrix(train_data[,3:ncol(train_data)])
y <- as.matrix(cbind(time=train_time,status=train_event))
t <- mean(train_data$time)

auc <- arc(filename, regtype="Lasso", metric="AUC", nfolds=5, standardize=T, alpha=.5, lambda=.01 )
print("Lasso AUC results")
print(auc)
auc <- arc(filename, regtype="Enet", metric="AUC", nfolds=5, alpha=.5, lambda=.01, standardize=T)
print("Enet AUC results")
print(auc)
auc <- arc(filename, regtype="Lapnet", metric="AUC", nfolds=5, alpha=.5, lambda=.1, standardize=T )
print("Lapnet AUC results")
print(auc)
auc <- arc(filename, regtype="Kernelnet", metric="AUC", nfolds=5, alpha=.5, lambda=.01, standardize=T )
print("Kernelnet AUC results")
print(auc)

mse <- arc(filename, regtype="Lasso", metric="MSE", nfolds=5, alpha=.5, lambda=.01, standardize=F )
print("Lasso MSE results")
print(mse)
mse <- arc(filename, regtype="Enet", metric="MSE", nfolds=5, alpha=.5, lambda=.01, standardize=T )
print("Enet MSE results")
print(mse)
mse <- arc(filename, regtype="Lapnet", metric="MSE", nfolds=5, alpha=.5, lambda=.01, standardize=T )
print("Lapnet MSE results")
print(mse)
mse <- arc(filename, regtype="Kernelnet", metric="MSE", nfolds=5, alpha=.5, lambda=.01, standardize=T )
print("Kernelnet MSE results")
print(mse)
