#Read TCGA BRCA expression data - log2(FPKM+1)
brca <- read.delim("C:/Users/nath0050/Desktop/BRCA/TCGA-BRCA.htseq_fpkm.tsv")
r1 <- strtrim(as.character(brca[,1]), 15)
c1 <- colnames(brca[2:ncol(brca)])
brca_exp <- as.matrix(brca[, 2:ncol(brca)])
rownames(brca_exp) <- r1
colnames(brca_exp) <- c1
rm(brca, r1, c1)

save(brca_exp, file="C:/Users/nath0050/Desktop/BRCA/TCGA_BRCA.RData")

#Normal vs. tumor
bc_n <- grep(".11", colnames(brca_exp))
bc_t <- grep(".01", colnames(brca_exp))

#Reducing genes from CTRP data
load("C:/Users/nath0050/Desktop/lncRNA/lncRNA_Rx/ctrp.auc.RData")
load("C:/Users/nath0050/Desktop/lncRNA/lncRNA_Rx/ctrp.rna.RData")

library(cluster)
tx <- proc.time()
clara.ctrp <- clara(t(ctrp.rna), k=200, metric="euclidean", stand=TRUE, keep.data = FALSE)
proc.time() - tx #248.35 

#Run ridge 
library(glmnet)
grep("birinapant", colnames(ctrp.auc)) #330 = "navitoclax:birinapant (1:1 mol/mol)", 410 = "birinapant"

yi <- which(!is.na(ctrp.auc[,410]))
auc.birin <- ctrp.auc[yi,410]
xi <- match(names(auc.birin), rownames(ctrp.rna))
ALL.birin <- ctrp.rna[xi, ]

names(auc.birin) <- rem_punc(names(auc.birin))
rownames(ALL.birin) <- rem_punc(rownames(ALL.birin))


#all
X <- scale(ALL.birin)
cv.all <- cv.glmnet(x = X, auc.birin, alpha=0)
kcv.all <- kfoldCV(X, auc.birin, k=10, a=0, lambda=cv.all$lambda.1se)

#vars
ALL.vars <- apply(ALL.birin, 2, var)
X <- scale(ALL.birin[, which(ALL.vars > 1)])
cv.var <- cv.glmnet(x = X, auc.birin, alpha=0)
kcv.var <- kfoldCV(X, auc.birin, k=10, a=0, lambda=cv.var$lambda.1se)

#cors
birin.cors <- apply(ALL.birin, 2, function (x) cor(x, auc.birin, method="spearman"))
X <- scale(ALL.birin[, which(abs(birin.cors) > 0.2)])
cv.cors <- cv.glmnet(x = X, auc.birin, alpha=0)
kcv.cors <- kfoldCV(X, auc.birin, k=10, a=0, lambda=cv.cors$lambda.1se)

library(randomForest)
krf.cors <- kRFCV(X, auc.birin, k=10, ntree=100)

#medoids 
X <- scale(t(clara.ctrp$medoids[,xi]))
cv.medoids <- cv.glmnet(x = X, auc.birin, alpha=0)
kcv.medoids <- kfoldCV(X, auc.birin, k=10, a=0, lambda=cv.medoids$lambda.1se)

#lasso + rf 
X <- scale(ALL.birin)
cv.la <- cv.glmnet(x = X, auc.birin, alpha=1)

#svm
library(e1071)
X <- scale(ALL.birin[, which(abs(birin.cors) > 0.2)])
cv.svm.cors <- kSVMCV(X, auc.birin, k=10)
X <- scale(ALL.birin[, which(ALL.vars > 1)])
cv.svm.var <- kSVMCV(X, auc.birin, k=10)

##### CCLE Feature Set ####
x2 <- match(names(auc.birin), colnames(CCLE.features))
auc.birin2 <- auc.birin[which(!is.na(x2))]
ALL.birin2 <- ALL.birin[which(!is.na(x2)), ]
ALL.feat <- t(CCLE.features[, na.omit(x2)])

feat.cors <- apply(ALL.feat, 2, function (x) cor(x, auc.birin2, method="pearson"))
X <- ALL.feat[, which(abs(feat.cors) > 0.15)]
krf.featcors <- kRFCV(X, auc.birin2, k=10, ntree=100)

birin2.cors <- apply(ALL.birin2, 2, function (x) cor(x, auc.birin2, method="spearman"))
X <- scale(ALL.birin2[, which(abs(birin2.cors) > 0.2)])
krf.birin2cors <- kRFCV(X, auc.birin2, k=10, ntree=100)

X <- cbind(ALL.feat[, which(abs(feat.cors) > 0.15)], ALL.birin2[, which(abs(birin2.cors) > 0.2)])
krf.corfeat <- kRFCV(X, auc.birin2, k=10, ntree=100)


#### Cluster AUC values ####
library(ROCR)
k2 <- kmeans(auc.birin, 2)

#cors
X <- scale(ALL.birin[, which(abs(birin.cors) > 0.2)])
krf.k2 <- kRFCV(X, k2$cluster, k=10, ntree=100)
temp <- randomForest(X, as.factor(k2$cluster), importance=TRUE, proximity=TRUE)
predictions <- as.vector(temp$votes[,2])
pred <- prediction(predictions, as.factor(k2$cluster))

perf_AUC <- performance(pred,"auc") #Calculate the AUC value
AUC <- perf_AUC@y.values[[1]]

perf_ROC <- performance(pred,"tpr","fpr") #plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5, 0.5, paste("AUC = ", format(AUC, digits=5, scientific=FALSE)))

#var
X <- scale(ALL.birin[, which(ALL.vars > 1)])
temp <- randomForest(X, as.factor(k2$cluster), importance=TRUE, proximity=TRUE)
predictions <- as.vector(temp$votes[,2])
pred <- prediction(predictions, as.factor(k2$cluster))

perf_AUC <- performance(pred,"auc") #Calculate the AUC value
AUC <- perf_AUC@y.values[[1]]

perf_ROC <- performance(pred,"tpr","fpr") #plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5, 0.5, paste("AUC = ", format(AUC, digits=5, scientific=FALSE)))


#cors + feat
X <- cbind(ALL.feat[, which(abs(feat.cors) > 0.15)], ALL.birin2[, which(abs(birin2.cors) > 0.2)])
k22 <- kmeans(auc.birin2, 2)
temp <- randomForest(X, as.factor(k22$cluster), importance=TRUE, proximity=TRUE)
predictions <- as.vector(temp$votes[,2])
pred <- prediction(predictions, as.factor(k22$cluster))

perf_AUC <- performance(pred,"auc") #Calculate the AUC value
AUC <- perf_AUC@y.values[[1]]

perf_ROC <- performance(pred,"tpr","fpr") #plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5, 0.5, paste("AUC = ", format(AUC, digits=5, scientific=FALSE)))

