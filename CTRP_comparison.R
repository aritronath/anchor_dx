load("CTRP_CV_results.RData")

dim(clara.ctrp$medoids) #200

ctrp.cors <- sapply(ctrp.comb_IKres, '[', 2)
ctrp.pvals <- sapply(ctrp.comb_IKres, '[', 3)

#Good prediction drugs 
hi <- which(ctrp.cors > 0.635) #top 10
colnames(ctrp.auc)[hi] 

#Density plot of AUC of good prediction drugs
par(bty='n', cex=1.5, mfrow=c(5,2))
apply(ctrp.auc[, hi], 2, function (x) {
  plot(density(x, na.rm=TRUE), main="") 
  abline(v=median(x, na.rm=TRUE))
  })

library(randomForest)
library(glmnet)
library(e1071)
#1. Feature selection + linear regression 
#A. Modeling response with genes selected by correlation + ridge regression
results.1A <- list()
tx <- proc.time()
for (i in 1:length(hi)) {
  y.auc <-  na.omit(ctrp.auc[, hi[i]])
  x.gex <- ctrp.rna[match(names(y.auc), rownames(ctrp.rna)), ]
  filt.gex <- corf(x.gex, y.auc, 200)
  
  results.1A[[i]] <- kGLMNETCV(filt.gex, y.auc, k=5, a=0, method="pearson")
  print(i)
}
time.1A <- proc.time() - tx

#B. Modeling response with LASSO
results.1B <- list()
tx <- proc.time()
for (i in 1:length(hi)) {
  y.auc <-  na.omit(ctrp.auc[, hi[i]])
  x.gex <- ctrp.rna[match(names(y.auc), rownames(ctrp.rna)), ]
  filt.gex <- corf(x.gex, y.auc, 200)
  
  results.1B[[i]] <- kGLMNETCV(filt.gex, y.auc, k=5, a=1, method="pearson")
  print(i)
}
time.1B <- proc.time() - tx

#2. Feature selection + non-linear regression 
#A. Modeling response with genes selected by correlation + random forests
results.2A <- list()
tx <- proc.time()
for (i in 1:length(hi)) {
  y.auc <-  na.omit(ctrp.auc[, hi[i]])
  x.gex <- ctrp.rna[match(names(y.auc), rownames(ctrp.rna)), ]
  filt.gex <- corf(x.gex, y.auc, 200)
  
  results.2A[[i]] <- kRFCV(filt.gex, y.auc, k=5, method="pearson")
  print(i)
}
time.2A <- proc.time() - tx


#B. Modeling response with genes selected by correlation + svm-r
results.2B <- list()
tx <- proc.time()
for (i in 1:length(hi)) {
  y.auc <-  na.omit(ctrp.auc[, hi[i]])
  x.gex <- ctrp.rna[match(names(y.auc), rownames(ctrp.rna)), ]
  filt.gex <- corf(x.gex, y.auc, 200)
  
  results.2B[[i]] <- kSVMCV(filt.gex, y.auc, k=5, method="pearson")
  print(i)
}
time.2B <- proc.time() - tx

#3. Dimensiion reduction + ridge regression 
results.3 <- list()
ctrp.pca <- svd(ctrp.rna, nu=10, nv=0)$u
rownames(ctrp.pca) <- rownames(ctrp.rna)
tx <- proc.time()
for (i in 1:length(hi)) {
  y.auc <-  na.omit(ctrp.auc[, hi[i]])
  filt.gex <- ctrp.pca[match(names(y.auc), rownames(ctrp.pca)), ]
  
  results.3[[i]] <- kGLMNETCV(filt.gex, y.auc, k=5, a=0, method="pearson")
  print(i)
}
time.3 <- proc.time() - tx

#4. Modeling response with universal medoids 
x.gex <- t(clara.ctrp$medoids)
rownames(x.gex) <- rownames(ctrp.rna)

#A. ridge regression 
results.4A <- list()
tx <- proc.time()
for (i in 1:length(hi)) {
  y.auc <-  na.omit(ctrp.auc[, hi[i]])
  filt.gex <- x.gex[match(names(y.auc), rownames(x.gex)), ]
  
  results.4A[[i]] <- kGLMNETCV(filt.gex, y.auc, k=5, a=0, method="pearson")
  print(i)
}
time.4A <- proc.time() - tx

#B. lasso
results.4B <- list()
tx <- proc.time()
for (i in 1:length(hi)) {
  y.auc <-  na.omit(ctrp.auc[, hi[i]])
  filt.gex <- x.gex[match(names(y.auc), rownames(x.gex)), ]
  
  results.4B[[i]] <- kGLMNETCV(filt.gex, y.auc, k=5, a=1, method="pearson")
  print(i)
}
time.4B <- proc.time() - tx

#C. random forests
results.4C <- list()
tx <- proc.time()
for (i in 1:length(hi)) {
  y.auc <-  na.omit(ctrp.auc[, hi[i]])
  filt.gex <- x.gex[match(names(y.auc), rownames(x.gex)), ]
  
  results.4C[[i]] <- kRFCV(filt.gex, y.auc, k=5, method="pearson")
  print(i)
}
time.4C <- proc.time() - tx

#D. svm-r
results.4D <- list()
tx <- proc.time()
for (i in 1:length(hi)) {
  y.auc <-  na.omit(ctrp.auc[, hi[i]])
  filt.gex <- x.gex[match(names(y.auc), rownames(x.gex)), ]
  
  results.4D[[i]] <- kSVMCV(filt.gex, y.auc, k=5, method="pearson")
  print(i)
}
time.4D <- proc.time() - tx

save(results.1A, results.1B, results.2A, results.2B, results.3, results.4A, results.4B, results.4C, results.4D,
     file="CTRP.Comparison.Results.RData")

save(time.1A, time.1B, time.2A, time.2B, time.3, time.4A, time.4B, time.4C, time.4D,
     file="CTRP.Comparison.time.RData")
