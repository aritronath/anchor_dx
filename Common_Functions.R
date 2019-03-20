#### Function to remove punctuations, switch to all uppercase ####
rem_punc <- function (dat) {
  temp <- toupper(gsub("\\-", "", dat))
  temp <- gsub("\\.", "", temp)
  temp <- gsub("\\_", "", temp)
  temp <- gsub(" ", "", temp)
  return(temp)
}  

### cutoff % of all samples 
posexp <- function (dat,cutoff) {
  index=0
  cut=cutoff*nrow(dat)/100
  for (i in 1:ncol(dat)) {
    if (sum(dat[,i]>1)>cut) 
      index=append(index,i)
  }
  index=index[2:length(index)]
  dat.f <- dat[,index]
  return (dat.f)
}


### Function to return matched matrices with common rows and columns
match_mat <- function (mat1,mat2) {  
  x <- match (rownames(mat1),rownames(mat2))
  x1 <- which(!is.na(x))
  x2 <- na.omit(x)
  
  y <- match (colnames(mat1),colnames(mat2))
  y1 <- which(!is.na(y))
  y2 <- na.omit(y)
  
  new.mat1 <- mat1[x1,y1]
  new.mat2 <- mat2[x2,y2]
  cmat <- list(new.mat1,new.mat2)
  return(cmat)
}


##Correlation filter
corf <- function (X, Y, num=200, method="pearson") {
  pcor <- apply(X, 2, function (x) abs(cor(x, Y, method=method)))
  gin <- which(rank(-pcor) < (num+1))
  filt_X <- X[, gin]
  return(filt_X)
}

##Variance filter
var.filter <-  function (dat, cutoff=0) {
  vars <- apply(dat, 1, var)
  dat <- dat[-which(vars == cutoff), ]
  return(dat)
}


#### Function for k-fold glmnet CV; LASSO a=1; Ridge a=0 ####
#depends on library(glmnet)
kGLMNETCV <- function (datX, datY, k=5, a=1, method="spearman") { 
  cv.fit <- cv.glmnet(datX, datY, nfolds=3, alpha=a)
  
  kgrp <- split(1:length(datY), sample(1:k, length(datY), replace=T))
  rmse <- array(dim=k); scor <- array(dim=k); pcor <- array(dim=k)
  
  for (i in 1:k) {
    testi <- unlist(kgrp[[i]])
    
    fit1 <- glmnet(datX[-testi, ], datY[-testi], alpha=a, lambda=cv.fit$lambda.1se)
    pred1 <- predict.glmnet(fit1, datX[testi, ])
    
    rmse[i] <- sqrt(sum((pred1 - datY[testi])^2)/length(pred1))
    cort <- cor.test(pred1, datY[testi], method=method)
    scor[i] <- cort$estimate
    pcor[i] <- cort$p.value
  }
  mean_err <- mean(rmse)
  mean_scor <- mean(scor, na.rm=TRUE)
  mean_pcor <- mean(pcor, na.rm=TRUE)
  return(cbind(mean_err, mean_scor, mean_pcor))
}

#### Function to perform K-fold RandomForest CV ####
#depends on library(randomForest)
kRFCV <- function (datX, datY, k=5, ntree=100, method="spearman") { 
  kgrp <- split(1:length(datY), sample(1:k, length(datY), replace=T))
  rmse <- array(dim=k); scor <- array(dim=k); pcor <- array(dim=k)
  
  for (i in 1:k) {
    testi <- unlist(kgrp[[i]])
    
    fit1 <- randomForest(datX[-testi, ], datY[-testi], ntree=ntree)
    pred1 <- predict(fit1, datX[testi, ])
    
    rmse[i] <- sqrt(sum((pred1 - datY[testi])^2)/length(pred1))
    cort <- cor.test(pred1, datY[testi], method=method)
    scor[i] <- cort$estimate
    pcor[i] <- cort$p.value
  }
  mean_err <- mean(rmse)
  mean_scor <- mean(scor, na.rm=TRUE)
  mean_pcor <- mean(pcor, na.rm=TRUE)
  return(cbind(mean_err, mean_scor, mean_pcor))
}


#### Function to perform K-fold SVM CV ####
#depends on library(e1071)
kSVMCV <- function (datX, datY, k=5, method="spearman") { 
  kgrp <- split(1:length(datY), sample(1:k, length(datY), replace=T))
  rmse <- array(dim=k); scor <- array(dim=k); pcor <- array(dim=k)
  
  for (i in 1:k) {
    testi <- unlist(kgrp[[i]])
    
    fit1 <- svm(datX[-testi, ], datY[-testi])
    pred1 <- predict(fit1, datX[testi, ])
    
    rmse[i] <- sqrt(sum((pred1 - datY[testi])^2)/length(pred1))
    cort <- cor.test(pred1, datY[testi], method=method)
    scor[i] <- cort$estimate
    pcor[i] <- cort$p.value
  }
  mean_err <- mean(rmse)
  mean_scor <- mean(scor, na.rm=TRUE)
  mean_pcor <- mean(pcor, na.rm=TRUE)
  return(cbind(mean_err, mean_scor, mean_pcor))
}

