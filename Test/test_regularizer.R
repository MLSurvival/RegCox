install_github("dmkdcs/RegCox/RegCox")
library(RegCox)

# Set the path to the file below
train_data <- read.csv("synthetic1-survival_norm.csv")
#train_data <- read.csv("breast_norm.csv")

train_time <- train_data$time
train_event <-train_data$event
x <- as.matrix(train_data[,3:ncol(train_data)])
y <- as.matrix(cbind(time=train_time,status=train_event))
t <- mean(train_data$time)

fit <- regcox(x, y, t, "Lasso",nfolds=5, stdbeta=F, standardize=F)
print("Lasso Results")
print(fit)
fit <- regcox(x, y, t, "Enet",nfolds=5 , alpha=0.1, stdbeta=F, standardize=F)
print("Enet Results")
print(fit)
fit <- regcox(x, y, t, "Lapnet", nfolds=5, alpha=0.1, stdbeta=F, standardize=F)
print("Lapnet Results")
print(fit)
fit <- regcox(x, y, t, "Kernelnet", nfolds=5, alpha=0.1, stdbeta=F, standardize=F)
print("Kernelnet Results")
print(fit)
fit <- regcox(x, y, t, "Oscar", nfolds=5, lambda1=0.01, lambda2=0.01, stdbeta=T, standardize=T)
print("Oscar Results")
print(fit)
fit <- regcox(x, y, t, "Fear", nfolds=5, tol=1e-3, lambda=0.1, itermax=100, stdbeta=T, standardize=F)
print("Fear Results")
print(fit)
fit <- regcox(x, y, t, "Cocktail", nfolds=5, tol=1e-3, lambda=0.1, itermax=100,alpha=.5, weight=1, stdbeta=T, standardize=F)
print("Cocktail Results")
print(fit)

