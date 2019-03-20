fulldist <- dist(t(fullmat), method = "euclidean")
full_hclust <- hclust(fulldist, method="complete")
save(fulldist, full_hclust, file="fullmat_clust.RData")

plot(full_hclust, labels=FALSE, col="grey30")
abline(h=45, lwd=2, col="blue")
clust_45 <- cutree(full_hclust, h=45)
clust_50 <- cutree(full_hclust, h=50)

#load gdsc gene expression data
load("rna")

#make distance/similarity matrix  
library(parallelDist)
gdscdist <- parDist(t(scale(rna)), method="euclidean", threads=6)
ccledist <- dist(t(scale(ccle.rna)), method="euclidean")


#get hierarchical clusters from gdsc expression 
gdsc_hclust <- hclust(gdscdist, method="complete")
save(gdscdist, gdsc_hclust, file="gdsc_clust.RData")

#plot clusters with lines at tree cut points
plot(gdsc_hclust, axes=FALSE, labels=FALSE, col="grey30", main="")
axis(side=2, at=seq(42.5, 55, 2.5), col="black", lwd=2, las=1, cex.axis=0.75)
abline(h=42.5, lwd=1.5, col="lightgreen", lty=2)
abline(h=45, lwd=1.5, col="lightgreen", lty=2)
abline(h=47.5, lwd=1.5, col="lightgreen", lty=2)
abline(h=50, lwd=1.5, col="lightgreen", lty=2)
abline(h=52.5, lwd=1.5, col="lightgreen", lty=2)
abline(h=55, lwd=1.5, col="lightgreen", lty=2)

cl_sizes=c("N=7", "N=22", "N=90", "N=268", "N=863", "N=2106")

#get clusters by cutting tree at different points #number of clusters
clg_40 <- cutree(gdsc_hclust, h=40) #dim(table(clg_40)) 4315
clg_42.5 <- cutree(gdsc_hclust, h=42.5) #2106
clg_45 <- cutree(gdsc_hclust, h=45) #863
clg_47.5 <- cutree(gdsc_hclust, h=47.5) #268
clg_50 <- cutree(gdsc_hclust, h=50) #90
clg_52.5 <- cutree(gdsc_hclust, h=52.5) #22
clg_55 <- cutree(gdsc_hclust, h=55) #7
clg_60 <- cutree(gdsc_hclust,h=60) #2


#function to extract 1 PC from each cluster of expression matrix 
get_pc1 <- function (expdat, treecut) {
  numc <- length(unique(treecut))
  pcmat <- matrix(data=NA, nrow=nrow(expdat), ncol=numc)
  for (i in 1:numc) {
    y <- which(treecut==i)
    selmat <- expdat[,y]
    pc1 <- svd(selmat, nu=1, nv=0)$u
    pcmat[,i] <- pc1
  }
  return(pcmat)
}

#Create matrix of PCs for each tree cut 
pc55 <- get_pc1(rna, clg_55)
pc52.5 <- get_pc1(rna, clg_52.5)
pc50 <- get_pc1(rna, clg_50)
pc47.5 <- get_pc1(rna, clg_47.5)
pc45 <- get_pc1(rna, clg_45)
pc42.5 <- get_pc1(rna, clg_42.5)

#function to obtain accuracy of predicting AUC
cvres <- function (cv.fit) {
  bf <- which(cv.fit$lambda.min == cv.fit$lambda)
  return (list(cv.fit$cvm[bf], cv.fit$cvsd[bf], cv.fit$nzero[bf]))
}

library(glmnet)
library(foreach)
library(doParallel)

#### Run glmnet ridge on parallel cores ####
cores=detectCores()
cl <- makeCluster(cores[1]-2)
registerDoParallel(cl)

#all genes
gdsc_ridge <- list(); gdsc_res <- list();
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match (names(auc), rownames(rna))
  rexp <- rna[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=0, nfolds=10, family="gaussian", standardize.response=TRUE, 
                       intercept=FALSE, parallel=TRUE)
  gdsc_ridge[[i]] <- coef (glmnet(rexp, auc, alpha=0, family="gaussian", standardize.response=TRUE, 
                                  intercept=FALSE, lambda=cv.fit1$lambda.min))
  gdsc_res[[i]] <- cvres(cv.fit1)
  
  print(i)
}  
save(gdsc_ridge, gdsc_res, file="GDSC_Ridge_Results.RData")

#pc55
pc55_ridge <- list(); pc55_res <- list()
rownames(pc55) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc55))
  rexp <- pc55[x,]

  cv.fit1 <- cv.glmnet(rexp, auc, alpha=0, nfolds=10, standardize.response=TRUE, 
                       intercept=FALSE, parallel=TRUE)
  pc55_ridge[[i]] <- coef (glmnet(rexp, auc, alpha=0, standardize.response=TRUE, 
                                  intercept=FALSE, lambda=cv.fit1$lambda.min))
  pc55_res[[i]] <- cvres(cv.fit1)
  
  print (i)  
}  

#pc52.5
pc52.5_ridge <- list(); pc52.5_res <- list()
rownames(pc52.5) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc52.5))
  rexp <- pc52.5[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=0, nfolds=10, standardize.response=TRUE,
                       intercept=FALSE, parallel=TRUE)
  pc52.5_ridge[[i]] <- coef (glmnet(rexp, auc, alpha=0, standardize.response=TRUE,
                                    intercept=FALSE, lambda=cv.fit1$lambda.min))
  pc52.5_res[[i]] <- cvres(cv.fit1)
  
  print (i)  
}  

#pc50
pc50_ridge <- list(); pc50_res <- list()
rownames(pc50) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc50))
  rexp <- pc50[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=0, nfolds=10, standardize.response=TRUE,
                       intercept=FALSE, parallel=TRUE)
  pc50_ridge[[i]] <- coef (glmnet(rexp, auc, alpha=0, standardize.response=TRUE,
                                  intercept=FALSE, lambda=cv.fit1$lambda.min))
  pc50_res[[i]] <- cvres(cv.fit1)
  
  print (i)  
}  

#pc47.5
pc47.5_ridge <- list(); pc47.5_res <- list()
rownames(pc47.5) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc47.5))
  rexp <- pc47.5[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=0, nfolds=10, standardize.response=TRUE,
                       intercept=FALSE, parallel=TRUE)
  pc47.5_ridge[[i]] <- coef (glmnet(rexp, auc, alpha=0, standardize.response=TRUE,
                                    intercept=FALSE, lambda=cv.fit1$lambda.min))
  pc47.5_res[[i]] <- cvres(cv.fit1)
  
  print (i)  
}  

#pc45
pc45_ridge <- list(); pc45_res <- list()
rownames(pc45) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc45))
  rexp <- pc45[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=0, nfolds=10, standardize.response=TRUE,
                       intercept=FALSE, parallel=TRUE)
  pc45_ridge[[i]] <- coef (glmnet(rexp, auc, alpha=0, standardize.response=TRUE,
                                  intercept=FALSE, lambda=cv.fit1$lambda.min))
  pc45_res[[i]] <- cvres(cv.fit1)
  
  print (i)  
}  

#pc42.5
pc42.5_ridge <- list(); pc42.5_res <- list()
rownames(pc42.5) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc42.5))
  rexp <- pc42.5[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=0, nfolds=10, standardize.response=TRUE,
                       intercept=FALSE, parallel=TRUE)
  pc42.5_ridge[[i]] <- coef (glmnet(rexp, auc, alpha=0, standardize.response=TRUE,
                                    intercept=FALSE, lambda=cv.fit1$lambda.min))
  pc42.5_res[[i]] <- cvres(cv.fit1)
  
  print (i)  
}  

stopCluster(cl)

save (pc42.5_ridge, pc42.5_res, pc45_ridge, pc45_res, pc47.5_ridge, pc47.5_res,
      pc50_ridge, pc50_res, pc52.5_ridge, pc52.5_res, pc55_ridge, pc55_res, 
      file="PC_Ridge_Results.RData")

pc42.5_m <- unlist(sapply(pc42.5_res, '[', 1))
pc45_m <- unlist(sapply(pc45_res, '[', 1))
pc47.5_m <- unlist(sapply(pc47.5_res, '[', 1))
pc50_m <- unlist(sapply(pc50_res, '[', 1))
pc52.5_m <- unlist(sapply(pc52.5_res, '[', 1))
pc55_m <- unlist(sapply(pc55_res, '[', 1))
gdsc_m <- unlist(sapply(gdsc_res, '[', 1))

boxplot(gdsc_m, pc55_m, pc52.5_m, pc50_m, pc47.5_m, pc45_m, pc42.5_m)

plot(gdsc_m)
points(pc55_m, pch=3, col=rainbow(1))
points(pc42.5_m, pch=4, col=rainbow(2))

t.test(gdsc_m, pc55_m)
t.test(gdsc_m, pc42.5_m)


#### Similar as above, custom CV function for COR and RMSE ####
#all genes
gdsc_Kres <- list();
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match (names(auc), rownames(rna))
  rexp <- rna[x,]

  gdsc_Kres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=FALSE, sr=TRUE)
  
  print(i)
}  

#pc55
pc55_Kres <- list()
rownames(pc55) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc55))
  rexp <- pc55[x,]
  
  pc55_Kres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=FALSE, sr=TRUE)
  
  print (i)  
}  

#pc52.5
pc52.5_Kres <- list()
rownames(pc52.5) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc52.5))
  rexp <- pc52.5[x,]
  
  pc52.5_Kres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=FALSE, sr=TRUE)

  print (i)  
}  

#pc50
pc50_Kres <- list()
rownames(pc50) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc50))
  rexp <- pc50[x,]
  
  pc50_Kres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=FALSE, sr=TRUE)
  
  print (i)  
}  

#pc47.5
pc47.5_Kres <- list()
rownames(pc47.5) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc47.5))
  rexp <- pc47.5[x,]

  pc47.5_Kres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=FALSE, sr=TRUE)
  
  print (i)  
}  

#pc45
pc45_Kres <- list()
rownames(pc45) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc45))
  rexp <- pc45[x,]
  
  pc45_Kres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=FALSE, sr=TRUE)
  
  print (i)  
}  

#pc42.5
pc42.5_Kres <- list()
rownames(pc42.5) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc42.5))
  rexp <- pc42.5[x,]
  
  pc42.5_Kres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=FALSE, sr=TRUE)
  
  print (i)  
}

save(gdsc_Kres, pc55_Kres, pc52.5_Kres, pc50_Kres, pc47.5_Kres, pc45_Kres, pc42.5_Kres, file="Custom_CV_gdsc_PC.RData")

#Unlist and plot stuff
pc55_Kr <- unlist(sapply(pc55_Kres, '[', 3))
pc52.5_Kr <- unlist(sapply(pc52.5_Kres, '[', 3))
pc50_Kr <- unlist(sapply(pc50_Kres, '[', 3))
pc47.5_Kr <- unlist(sapply(pc47.5_Kres, '[', 3))
pc45_Kr <- unlist(sapply(pc45_Kres, '[', 3))
pc42.5_Kr <- unlist(sapply(pc42.5_Kres, '[', 3))
gdsc_Kr <- unlist(sapply(gdsc_Kres, '[', 3))

pc55_Kg <- unlist(sapply(pc55_Kres, '[', 1))
pc52.5_Kg <- unlist(sapply(pc52.5_Kres, '[', 1))
pc50_Kg <- unlist(sapply(pc50_Kres, '[', 1))
pc47.5_Kg <- unlist(sapply(pc47.5_Kres, '[', 1))
pc45_Kg <- unlist(sapply(pc45_Kres, '[', 1))
pc42.5_Kg <- unlist(sapply(pc42.5_Kres, '[', 1))
gdsc_Kg <- unlist(sapply(gdsc_Kres, '[', 1))

pc55_Kc <- unlist(sapply(pc55_Kres, '[', 2))
pc52.5_Kc <- unlist(sapply(pc52.5_Kres, '[', 2))
pc50_Kc <- unlist(sapply(pc50_Kres, '[', 2))
pc47.5_Kc <- unlist(sapply(pc47.5_Kres, '[', 2))
pc45_Kc <- unlist(sapply(pc45_Kres, '[', 2))
pc42.5_Kc <- unlist(sapply(pc42.5_Kres, '[', 2))
gdsc_Kc <- unlist(sapply(gdsc_Kres, '[', 2))

boxplot(pc55_Kr, pc52.5_Kr, pc50_Kr, pc47.5_Kr, pc45_Kr, pc42.5_Kr, gdsc_Kr, ylim=c(-1,1))
boxplot(pc55_Kg, pc52.5_Kg, pc50_Kg, pc47.5_Kg, pc45_Kg, pc42.5_Kg, gdsc_Kg)
boxplot(pc55_Kc, pc52.5_Kc, pc50_Kc, pc47.5_Kc, pc45_Kc, pc42.5_Kc, gdsc_Kc)


#### With Intercept, similar as above, custom CV function for COR and RMSE ####
#all genes
gdsc_IKres <- list();
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match (names(auc), rownames(rna))
  rexp <- rna[x,]
  
  gdsc_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print(i)
}  

#pc55
pc55_IKres <- list()
rownames(pc55) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc55))
  rexp <- pc55[x,]
  
  pc55_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}  

#pc52.5
pc52.5_IKres <- list()
rownames(pc52.5) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc52.5))
  rexp <- pc52.5[x,]
  
  pc52.5_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}  

#pc50
pc50_IKres <- list()
rownames(pc50) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc50))
  rexp <- pc50[x,]
  
  pc50_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}  

#pc47.5
pc47.5_IKres <- list()
rownames(pc47.5) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc47.5))
  rexp <- pc47.5[x,]
  
  pc47.5_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}  

#pc45
pc45_IKres <- list()
rownames(pc45) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc45))
  rexp <- pc45[x,]
  
  pc45_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}  

#pc42.5
pc42.5_IKres <- list()
rownames(pc42.5) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc42.5))
  rexp <- pc42.5[x,]
  
  pc42.5_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}

save(gdsc_IKres, pc55_IKres, pc52.5_IKres, pc50_IKres, pc47.5_IKres, pc45_IKres, pc42.5_IKres, file="Custom_Intercept_CV_gdsc_PC.RData")

#Unlist and plot stuff
pc55_IKr <- unlist(sapply(pc55_IKres, '[', 3))
pc52.5_IKr <- unlist(sapply(pc52.5_IKres, '[', 3))
pc50_IKr <- unlist(sapply(pc50_IKres, '[', 3))
pc47.5_IKr <- unlist(sapply(pc47.5_IKres, '[', 3))
pc45_IKr <- unlist(sapply(pc45_IKres, '[', 3))
pc42.5_IKr <- unlist(sapply(pc42.5_IKres, '[', 3))
gdsc_IKr <- unlist(sapply(gdsc_IKres, '[', 3))

pc55_IKg <- unlist(sapply(pc55_IKres, '[', 1))
pc52.5_IKg <- unlist(sapply(pc52.5_IKres, '[', 1))
pc50_IKg <- unlist(sapply(pc50_IKres, '[', 1))
pc47.5_IKg <- unlist(sapply(pc47.5_IKres, '[', 1))
pc45_IKg <- unlist(sapply(pc45_IKres, '[', 1))
pc42.5_IKg <- unlist(sapply(pc42.5_IKres, '[', 1))
gdsc_IKg <- unlist(sapply(gdsc_IKres, '[', 1))

pc55_IKc <- unlist(sapply(pc55_IKres, '[', 2))
pc52.5_IKc <- unlist(sapply(pc52.5_IKres, '[', 2))
pc50_IKc <- unlist(sapply(pc50_IKres, '[', 2))
pc47.5_IKc <- unlist(sapply(pc47.5_IKres, '[', 2))
pc45_IKc <- unlist(sapply(pc45_IKres, '[', 2))
pc42.5_IKc <- unlist(sapply(pc42.5_IKres, '[', 2))
gdsc_IKc <- unlist(sapply(gdsc_IKres, '[', 2))

boxplot(pc55_IKr, pc52.5_IKr, pc50_IKr, pc47.5_IKr, pc45_IKr, pc42.5_IKr, gdsc_IKr, CFE_IKr, iPC10_IKr)
boxplot(pc55_IKg, pc52.5_IKg, pc50_IKg, pc47.5_IKg, pc45_IKg, pc42.5_IKg, gdsc_IKg)
boxplot(pc55_IKc, pc52.5_IKc, pc50_IKc, pc47.5_IKc, pc45_IKc, pc42.5_IKc, gdsc_IKc)

ks.test(pc55_IKr, pc52.5_IKr)
ks.test(pc52.5_IKr, pc50_IKr)
ks.test(pc50_IKr, pc47.5_IKr)
ks.test(pc47.5_IKr, pc45_IKr)
ks.test(pc45_IKr, pc42.5_IKr)

ks.test(pc42.5_IKr, gdsc_IKr)
ks.test(pc45_IKr, gdsc_IKr)
ks.test(pc47.5_IKr, gdsc_IKr)

#### With Intercept and CFEs, similar as above, custom CV function for COR and RMSE ####
#This rearranges cell functional events according to 
#rauc - rownames in auc.gdsc, inherited from lncRNA_Rx_Ridge.R
x2 <- match(rauc, toupper(rownames(CFE)))
auc.CFE <- CFE[x2,]

#all genes
CFE_gdsc_IKres <- list();
for (i in 1:ncol(auc.gdsc)) {
  ax <- which(!is.na(auc.gdsc[,i]))
  auc <- auc.gdsc[ax, i]
  x <- match (names(auc), rownames(rna))
  rexp <- cbind(rna[x,], auc.CFE[ax,]) 
  
  CFE_gdsc_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print(i)
}  

# CFE by itself
CFE_IKres <- list()
for (i in 1:ncol(auc.gdsc)) {
  ax <- which(!is.na(auc.gdsc[,i]))
  auc <- auc.gdsc[ax, i]
  
  CFE_IKres[[i]] <- kfoldCV(auc.CFE[ax,], auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}

save(CFE_IKres, CFE_gdsc_IKres, file="Custom_Intercept_CV_gdsc_PC_CFE.RData")

CFE_gdsc_IKr <- unlist(sapply(CFE_gdsc_IKres, '[', 3))
CFE_IKr <- unlist(sapply(CFE_IKres, '[', 3))

boxplot(CFE_gdsc_IKr, CFE_IKr)
boxplot(gdsc_IKr, CFE_gdsc_IKr)

CFE_gdsc_IKr <- unlist(sapply(CFE_gdsc_IKres, '[', 3))
CFE_IKr <- unlist(sapply(CFE_IKres, '[', 3))

boxplot(CFE_gdsc_IKr, CFE_IKr)
boxplot(gdsc_IKr, CFE_gdsc_IKr)

#pc55 + CFE + TT
pc55_CT_IKres <- list()
for (i in 1:ncol(auc.gdsc)) {
  ax <- which(!is.na(auc.gdsc[,i]))
  auc <- auc.gdsc[ax, i]
  x <- match(names(auc), rownames(pc55))
  rexp <- cbind(pc55[x,], auc.CFE[ax,], as.numeric(auc.tiss[ax]))
  
  pc55_CT_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}  
pc55_CT_IKr <- unlist(sapply(pc55_CT_IKres, '[', 3))


#pc52.5 + CFE + TT
pc52.5_CT_IKres <- list()
for (i in 1:ncol(auc.gdsc)) {
  ax <- which(!is.na(auc.gdsc[,i]))
  auc <- auc.gdsc[ax, i]
  x <- match(names(auc), rownames(pc52.5))
  rexp <- cbind(pc52.5[x,], auc.CFE[ax,], as.numeric(auc.tiss[ax]))
  
  pc52.5_CT_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}  
pc52.5_CT_IKr <- unlist(sapply(pc52.5_CT_IKres, '[', 3))


#pc50 + CFE + TT
pc50_CT_IKres <- list()
for (i in 1:ncol(auc.gdsc)) {
  ax <- which(!is.na(auc.gdsc[,i]))
  auc <- auc.gdsc[ax, i]
  x <- match(names(auc), rownames(pc50))
  rexp <- cbind(pc50[x,], auc.CFE[ax,], as.numeric(auc.tiss[ax]))
  
  pc50_CT_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}  
pc50_CT_IKr <- unlist(sapply(pc50_CT_IKres, '[', 3))

boxplot(pc55_CT_IKr, pc55_IKr, pc52.5_CT_IKr, pc52.5_IKr, pc50_CT_IKr, pc50_IKr)


### 
pc50_T_IKres <- list()
for (i in 1:ncol(auc.gdsc)) {
  ax <- which(!is.na(auc.gdsc[,i]))
  auc <- auc.gdsc[ax, i]
  x <- match(names(auc), rownames(pc50))
  rexp <- cbind(pc50[x,], as.numeric(auc.tiss[ax]))
  
  pc50_T_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}  
pc50_T_IKr <- unlist(sapply(pc50_T_IKres, '[', 3))


#### Relationship with GDSC drugs: targeted vs. non-targeted ####
drugs.gdsc <- read.csv("C:/Users/nath0050/Desktop/GDSC/Drugs_GDSC.csv")
tc <- drugs.gdsc[,6]
boxplot(CFE_pc55_IKr~tc)
boxplot(gdsc_IKr~tc)

#### Number of samples/drug ####
gdsc.nsd <- array(dim=ncol(auc.gdsc))
for (i in 1:ncol(auc.gdsc)) {
  gdsc.nsd[i] <- sum(!is.na(auc.gdsc[,i]))
}

#### How different are clusters: Jaccard index ####
library(cluseval)
cluster_similarity(clg_40, clg_42.5, similarity="jaccard", method="independence")
cluster_similarity(clg_42.5, clg_45, similarity="jaccard", method="independence")
cluster_similarity(clg_45, clg_47.5, similarity="jaccard", method="independence")
cluster_similarity(clg_47.5, clg_50, similarity="jaccard", method="independence")
cluster_similarity(clg_50, clg_52.5, similarity="jaccard", method="independence")
cluster_similarity(clg_52.5, clg_55, similarity="jaccard", method="independence")

cluster_similarity(clg_40, clg_55, similarity="jaccard", method="independence")
cluster_similarity(clg_40, clg_52.5, similarity="jaccard", method="independence")
cluster_similarity(clg_40, clg_50, similarity="jaccard", method="independence")
cluster_similarity(clg_40, clg_47.5, similarity="jaccard", method="independence")
cluster_similarity(clg_40, clg_45, similarity="jaccard", method="independence")
cluster_similarity(clg_40, clg_42.5, similarity="jaccard", method="independence")

#### Get Ridge coefficients for clusters #### 
#PC55 N=7
pc55_ridge_coef <- matrix(nrow=ncol(auc.gdsc), ncol=ncol(pc55))
for (i in 1:ncol(auc.gdsc)) {
  ax <- which(!is.na(auc.gdsc[,i]))
  auc <- auc.gdsc[ax, i]
  x <- match (names(auc), rownames(pc55))
  rexp <- pc55[x,] 
  
  cv.fit <- cv.glmnet(rexp, auc, nfolds=5)
  ridge.fit <- glmnet(rexp, auc, alpha=0, lambda=cv.fit$lambda.1se)
  
  pc55_ridge_coef[i,] <- as.array(ridge.fit$beta)
  print(i)
} 
write.csv(pc55_ridge_coef, file="C:/Users/nath0050/Desktop/Reducing Genes/pc55_ridge_coef.csv")

#PC52.5 N=22
pc52.5_ridge_coef <- matrix(nrow=ncol(auc.gdsc), ncol=ncol(pc52.5))
for (i in 1:ncol(auc.gdsc)) {
  ax <- which(!is.na(auc.gdsc[,i]))
  auc <- auc.gdsc[ax, i]
  x <- match (names(auc), rownames(pc52.5))
  rexp <- pc52.5[x,] 
  
  cv.fit <- cv.glmnet(rexp, auc, nfolds=5)
  ridge.fit <- glmnet(rexp, auc, alpha=0, lambda=cv.fit$lambda.1se)
  
  pc52.5_ridge_coef[i,] <- as.array(ridge.fit$beta)
  print(i)
} 
write.csv(pc52.5_ridge_coef, file="C:/Users/nath0050/Desktop/Reducing Genes/pc52.5_ridge_coef.csv")

#PC50 N=90
pc50_ridge_coef <- matrix(nrow=ncol(auc.gdsc), ncol=ncol(pc50))
for (i in 1:ncol(auc.gdsc)) {
  ax <- which(!is.na(auc.gdsc[,i]))
  auc <- auc.gdsc[ax, i]
  x <- match (names(auc), rownames(pc50))
  rexp <- pc50[x,] 
  
  cv.fit <- cv.glmnet(rexp, auc, nfolds=5)
  ridge.fit <- glmnet(rexp, auc, alpha=0, lambda=cv.fit$lambda.1se)
  
  pc50_ridge_coef[i,] <- as.array(ridge.fit$beta)
  print(i)
} 
write.csv(pc50_ridge_coef, file="C:/Users/nath0050/Desktop/Reducing Genes/pc50_ridge_coef.csv")

#PC47.5 N=268
pc47.5_ridge_coef <- matrix(nrow=ncol(auc.gdsc), ncol=ncol(pc47.5))
for (i in 1:ncol(auc.gdsc)) {
  ax <- which(!is.na(auc.gdsc[,i]))
  auc <- auc.gdsc[ax, i]
  x <- match (names(auc), rownames(pc47.5))
  rexp <- pc47.5[x,] 
  
  cv.fit <- cv.glmnet(rexp, auc, nfolds=5)
  ridge.fit <- glmnet(rexp, auc, alpha=0, lambda=cv.fit$lambda.1se)
  
  pc47.5_ridge_coef[i,] <- as.array(ridge.fit$beta)
  print(i)
} 
write.csv(pc47.5_ridge_coef, file="C:/Users/nath0050/Desktop/Reducing Genes/pc47.5_ridge_coef.csv")

#
su_47.5 <- apply(pc47.5_ridge_coef, 2, median)
plot(su_47.5, cex=0.01)
text(labels=c(1:length(su_47.5)), x=c(1:length(su_47.5)), y=su_47.5)

su_50 <- apply(pc50_ridge_coef, 2, median)
plot(su_50, cex=0.01)
text(labels=c(1:length(su_50)), x=c(1:length(su_50)), y=su_50)

su_52.5 <- apply(pc52.5_ridge_coef, 2, median)
plot(su_52.5, cex=0.01)
text(labels=c(1:length(su_52.5)), x=c(1:length(su_52.5)), y=su_52.5)

su_55 <- apply(pc55_ridge_coef, 2, median)
plot(su_55, cex=0.01)
text(labels=c(1:length(su_55)), x=c(1:length(su_55)), y=su_55)

boxplot(pc50_ridge_coef)
boxplot(pc55_ridge_coef)

#### Cluster pathway enrichment analyses with cluster profiler ######
library(clusterProfiler)
library(org.Hs.eg.db)

#make list of modules with genes corresponding to those modules
temp_genes <- bitr(colnames(auc.gdsc.rna), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
x <- which(duplicated(temp_genes[,2]))
temp_genes <- temp_genes[-x,]

## For clg_47.5 (90 clusters) 
genesbymod_47.5 <- list()
cln <- unique(clg_47.5)
for (i in 1:length(cln)) {
  y <- which(clg_47.5 == i)
  x <- match(colnames(auc.gdsc.rna)[y], temp_genes[,1])
  genesbymod_47.5[[i]] <- temp_genes[x,2]
}
names(genesbymod_47.5) <- cln

ck_47.5 <- compareCluster(genesbymod_47.5, fun="enrichKEGG", organism = "hsa", keyType = "kegg", 
                     pvalueCutoff = 0.05, pAdjustMethod = "BH")
dotplot(ck_47.5)
ck.goall_47.5 <- compareCluster(genesbymod_47.5, fun="enrichGO", OrgDb = "org.Hs.eg.db", ont="ALL", 
                              pvalueCutoff = 0.01, pAdjustMethod = "BH", pool=TRUE)

## For clg_50 (90 clusters) 
genesbymod_50 <- list()
cln <- unique(clg_50)
for (i in 1:length(cln)) {
  y <- which(clg_50 == i)
  x <- match(colnames(auc.gdsc.rna)[y], temp_genes[,1])
  genesbymod_50[[i]] <- temp_genes[x,2]
}
names(genesbymod_50) <- cln

ck <- compareCluster(genesbymod_50, fun="enrichKEGG", organism = "hsa", keyType = "kegg", 
                     pvalueCutoff = 0.05, pAdjustMethod = "BH")
dotplot(ck)
ck.goall_50 <- compareCluster(genesbymod_50, fun="enrichGO", OrgDb = "org.Hs.eg.db", ont="ALL", 
                           pvalueCutoff = 0.01, pAdjustMethod = "BH", pool=TRUE)
dotplot(ck.goall)
save(genesbymod_50, ck, ck.goall, file="C:/Users/nath0050/Desktop/Reducing Genes/CLG50_ClusterProfiles.RData")

ck.mat <- fortify(ck)
ck.goall.mat <- fortify(ck.goall) 

write.csv(ck.mat, file="C:/Users/nath0050/Desktop/Reducing Genes/CLG50_KEGG.csv")
write.csv(ck.goall.mat, file="C:/Users/nath0050/Desktop/Reducing Genes/CLG50_GOALL.csv")


## For clg_52.5 (22 clusters) 
genesbymod_52.5 <- list()
cln <- unique(clg_52.5)
for (i in 1:length(cln)) {
  y <- which(clg_52.5 == i)
  x <- match(colnames(auc.gdsc.rna)[y], temp_genes[,1])
  genesbymod_52.5[[i]] <- temp_genes[x,2]
}
names(genesbymod_52.5) <- cln

ck_52.5 <- compareCluster(genesbymod_52.5, fun="enrichKEGG", organism = "hsa", keyType = "kegg", 
                     pvalueCutoff = 0.05, pAdjustMethod = "BH")
dotplot(ck_52.5)
ck.goall_52.5 <- compareCluster(genesbymod_52.5, fun="enrichGO", OrgDb = "org.Hs.eg.db", ont="ALL", 
                           pvalueCutoff = 0.01, pAdjustMethod = "BH", pool=TRUE)
dotplot(ck.goall_52.5)
save(genesbymod_52.5, ck_52.5, ck.goall_52.5, file="C:/Users/nath0050/Desktop/Reducing Genes/CLG52.5_ClusterProfiles.RData")

ck.mat_52.5 <- fortify(ck_52.5)
ck.goall.mat_52.5 <- fortify(ck.goall_52.5) 

write.csv(ck.mat_52.5, file="C:/Users/nath0050/Desktop/Reducing Genes/CLG52.5_KEGG.csv")
write.csv(ck.goall.mat_52.5, file="C:/Users/nath0050/Desktop/Reducing Genes/CLG52.5_GOALL.csv")


## For clg_55 (7 clusters) 
genesbymod_55 <- list()
cln <- unique(clg_55)
for (i in 1:length(cln)) {
  y <- which(clg_55 == i)
  x <- match(colnames(auc.gdsc.rna)[y], temp_genes[,1])
  genesbymod_55[[i]] <- temp_genes[x,2]
}
names(genesbymod_55) <- cln

ck_55 <- compareCluster(genesbymod_55, fun="enrichKEGG", organism = "hsa", keyType = "kegg", 
                          pvalueCutoff = 0.05, pAdjustMethod = "BH")
dotplot(ck_55)
ck.goall_55 <- compareCluster(genesbymod_55, fun="enrichGO", OrgDb = "org.Hs.eg.db", ont="ALL", 
                           pvalueCutoff = 0.01, pAdjustMethod = "BH", pool=TRUE)
dotplot(ck.goall_55)
save(genesbymod_55, ck_55, ck.goall_55, file="C:/Users/nath0050/Desktop/Reducing Genes/CLG55_ClusterProfiles.RData")

ck.mat_55 <- fortify(ck_55)
ck.goall.mat_55 <- fortify(ck.goall_55) 

write.csv(ck.mat_55, file="C:/Users/nath0050/Desktop/Reducing Genes/CLG55_KEGG.csv")
write.csv(ck.goall.mat_55, file="C:/Users/nath0050/Desktop/Reducing Genes/CLG55_GOALL.csv")
