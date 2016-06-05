#' @export

arc.lapnet.auc <- function(train_data1, nfolds, alpha, lambda)
{
  
  folds <- coxsplit(train_data1, nfolds)
  auc <- numeric(nfolds)
  all_results <- NULL
  
  for (i in 1:nfolds)
  {
    
    train_data <- train_data1[which(folds!=i),]
    test_data <- train_data1[which(folds==i),]
    
    size <- nrow(train_data)
    nfeatures <- ncol(train_data)
    
    init_train_size <- round(.1*nrow(train_data))
    init_size <- round(init_train_size/2)
    increment_size <- round(.2*init_train_size)
    iterations <- round((.6*nrow(train_data)-init_train_size)/increment_size)
    
    train_data_cens <- subset(train_data, train_data$event == 0)
    train_data_uncens <- subset(train_data, train_data$event == 1)
    
    size_cens <- nrow(train_data_cens)
    size_uncens <- nrow(train_data_uncens)
    
    train_time <- train_data_cens$time
    train_event <-train_data_cens$event
    x_cens <- as.matrix(train_data_cens[,3:nfeatures])
    y_cens <- as.matrix(cbind(time=train_time,status=train_event))
    
    train_time <- train_data_uncens$time
    train_event <-train_data_uncens$event
    x_uncens <- as.matrix(train_data_uncens[,3:nfeatures])
    y_uncens <- as.matrix(cbind(time=train_time,status=train_event))
    
    xtrain <- rbind(x_cens[1:init_size,], x_uncens[1:init_size,])
    xtest <- rbind(x_cens[(init_size+1):size_cens,], x_uncens[(init_size+1):size_uncens,])
    ytrain <- rbind(y_cens[1:init_size,], y_uncens[1:init_size,])
    ytest <- rbind(y_cens[(init_size+1):size_cens,], y_uncens[(init_size+1):size_uncens,])
    
    results <- data.frame(train_size=0, aucfold=0)
	results <- NULL
    
    t <- mean(train_data$time)
    
    for(j in c(1:iterations))
    {
      
	  Laplacian <- getLaplacian(xtrain)
  
	  fit <- Coxnet(xtrain, ytrain, Omega = Laplacian, penalty="Net", alpha=alpha, lambda=lambda, isd=TRUE)
	  beta <- fit$Beta
      
      dt <- data.frame(cbind(ytrain,xtrain))
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
      
      # find the predictions on the training data
      ht <- -log(s0^(exp(xtest%*%as.matrix(beta))))
      
      output <- concordance.index(x=ht, surv.time=ytest[,1], surv.event=ytest[,2], method="noether")
      aucfold <- output$c.index
      
	if(is.null(results))
    {
      results <- data.frame(train_size=dim(xtrain)[1], aucfold=aucfold)
    }
    else
    {
      results <- rbind(results, c(dim(xtrain)[1], aucfold))
    }
      
      ################## Find Least Predictions ############################
      
      ind <- sort(ht, index.return=T, decreasing=F)
      
      # Find least of sample size
      ind <- ind$ix[1:increment_size]
      xtrain <- rbind(xtrain, xtest[ind,])
      ytrain <- rbind(ytrain, ytest[ind,])
      xtest <- xtest[-ind,]
      ytest <- ytest[-ind,]
      
    }
 
    # find the prediction on the test data
    
    test_time <- test_data$time
    test_event <- test_data$event
    xtest <- as.matrix(test_data[,3:nfeatures])
    ytest <- as.matrix(cbind(time=test_time,status=test_event))
    
    ht <- -log(s0^(exp(xtest%*%as.matrix(beta))))
    
    output <- concordance.index(x=ht, surv.time=ytest[,1], surv.event=ytest[,2], method="noether")
    auc[i] <- output$c.index
    
    if(is.null(all_results))
    {
      all_results <- results
    }
    else
    {
      all_results <- cbind(all_results, results$aucfold)
    }
    
  }
  
  all_results <- cbind(train_size=all_results[,1], auc=rowSums(all_results[,c(2:(nfolds+1))])/nfolds)
  
  return(c(list(nfold_auc=auc), auc_mean = mean(auc), auc_sd=sd(auc), list(auc_seq=all_results)))
  
}