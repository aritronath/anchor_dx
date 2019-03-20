gdsc_lasso <- list(); gdsc_lasso_res <- list();
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match (names(auc), rownames(rna))
  rexp <- rna[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=1, nlambda=20, nfolds=10, family="gaussian", standardize.response=TRUE, 
                       intercept=FALSE, parallel=TRUE)
  gdsc_lasso[[i]] <- coef (glmnet(rexp, auc, alpha=1, family="gaussian", standardize.response=TRUE, 
                                  intercept=FALSE, lambda=cv.fit1$lambda.min))
  gdsc_lasso_res[[i]] <- cvres(cv.fit1)
  
  print(i)
}

save(gdsc_lasso, gdsc_lasso_res, file="GDSC_Lasso_Results.RData")

#pc55
pc55_lasso <- list(); pc55_lasso_res <- list()
rownames(pc55) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc55))
  rexp <- pc55[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=1, nlambda=20, nfolds=10, standardize.response=TRUE, 
                       intercept=FALSE, parallel=TRUE)
  pc55_lasso[[i]] <- coef (glmnet(rexp, auc, alpha=1, standardize.response=TRUE, 
                                  intercept=FALSE, lambda=cv.fit1$lambda.min))
  pc55_lasso_res[[i]] <- cvres(cv.fit1)
  
  print (i)  
}  

#pc52.5
pc52.5_lasso <- list(); pc52.5_lasso_res <- list()
rownames(pc52.5) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc52.5))
  rexp <- pc52.5[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=1, nlambda=20, nfolds=10, standardize.response=TRUE,
                       intercept=FALSE, parallel=TRUE)
  pc52.5_lasso[[i]] <- coef (glmnet(rexp, auc, alpha=1, standardize.response=TRUE,
                                    intercept=FALSE, lambda=cv.fit1$lambda.min))
  pc52.5_lasso_res[[i]] <- cvres(cv.fit1)
  
  print (i)  
}  

#pc50
pc50_lasso <- list(); pc50_lasso_res <- list()
rownames(pc50) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc50))
  rexp <- pc50[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=1, nfolds=10, nlambda=20, standardize.response=TRUE,
                       intercept=FALSE, parallel=TRUE)
  pc50_lasso[[i]] <- coef (glmnet(rexp, auc, alpha=1, standardize.response=TRUE,
                                  intercept=FALSE, lambda=cv.fit1$lambda.min))
  pc50_lasso_res[[i]] <- cvres(cv.fit1)
  
  print (i)  
}  

#pc47.5
pc47.5_lasso <- list(); pc47.5_lasso_res <- list()
rownames(pc47.5) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc47.5))
  rexp <- pc47.5[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=1, nfolds=10, standardize.response=TRUE,
                       intercept=FALSE, parallel=TRUE)
  pc47.5_lasso[[i]] <- coef (glmnet(rexp, auc, alpha=1, standardize.response=TRUE,
                                    intercept=FALSE, lambda=cv.fit1$lambda.min))
  pc47.5_lasso_res[[i]] <- cvres(cv.fit1)
  
  print (i)  
}  

#pc45
pc45_lasso <- list(); pc45_lasso_res <- list()
rownames(pc45) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc45))
  rexp <- pc45[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=1, nfolds=10, standardize.response=TRUE,
                       intercept=FALSE, parallel=TRUE)
  pc45_lasso[[i]] <- coef (glmnet(rexp, auc, alpha=1, standardize.response=TRUE,
                                  intercept=FALSE, lambda=cv.fit1$lambda.min))
  pc45_lasso_res[[i]] <- cvres(cv.fit1)
  
  print (i)  
}  

#pc42.5
pc42.5_lasso <- list(); pc42.5_lasso_res <- list()
rownames(pc42.5) <- rownames(rna)
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc42.5))
  rexp <- pc42.5[x,]
  
  cv.fit1 <- cv.glmnet(rexp, auc, alpha=1, nfolds=10, standardize.response=TRUE,
                       intercept=FALSE, parallel=TRUE)
  pc42.5_lasso[[i]] <- coef (glmnet(rexp, auc, alpha=1, standardize.response=TRUE,
                                    intercept=FALSE, lambda=cv.fit1$lambda.min))
  pc42.5_lasso_res[[i]] <- cvres(cv.fit1)
  
  print (i)  
}  

save (pc42.5_lasso, pc42.5_lasso_res, pc45_lasso, pc45_lasso_res, pc47.5_lasso, pc47.5_lasso_res,
      pc50_lasso, pc50_lasso_res, pc52.5_lasso, pc52.5_lasso_res, pc55_lasso, pc55_lasso_res, 
      file="PC_Lasso_Results.RData")


stopCluster(cl)

pc42.5_lm <- unlist(sapply(pc42.5_lasso_res, '[', 1))
pc45_lm <- unlist(sapply(pc45_lasso_res, '[', 1))
pc47.5_lm <- unlist(sapply(pc47.5_lasso_res, '[', 1))
pc50_lm <- unlist(sapply(pc50_lasso_res, '[', 1))
pc52.5_lm <- unlist(sapply(pc52.5_lasso_res, '[', 1))
pc55_lm <- unlist(sapply(pc55_lasso_res, '[', 1))
gdsc_lm <- unlist(sapply(gdsc_lasso_res, '[', 1))

boxplot(pc55_lm, pc52.5_lm, pc50_lm, pc47.5_lm, pc45_lm, pc42.5_lm, gdsc_lm)

