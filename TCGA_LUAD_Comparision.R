load("C:/Users/nath0050/Desktop/TCGA_lung/TCGA_lung.RData")
colnames(luad) <- gsub("\\.", "\\-", colnames(luad))

luad.pheno <- read.delim("C:/Users/nath0050/Desktop/TCGA_lung/TCGA-LUAD.GDC_phenotype.tsv", header=TRUE)
luad.surv <- read.delim("C:/Users/nath0050/Desktop/TCGA_lung/TCGA-LUAD.survival.tsv", header=TRUE)

#### Clean-up ####
var.filter <-  function (dat, cutoff=0) {
  vars <- apply(dat, 1, var)
  dat <- dat[-which(vars == cutoff), ]
  return(dat)
}
luad <- var.filter(luad)

luad <- luad[, -grep("11A", colnames(luad))]

#CLARA clustering for medoids
library(cluster)
x <- match(colnames(luad), luad.surv$sample)
luad.2 <- t(luad[,which(!is.na(x))])
luad.surv.2 <- luad.surv[na.omit(x),]

tx <- proc.time()
luad.clara <- clara(t(luad.2), k=200, metric="euclidean", stand=TRUE, keep.data=FALSE)
proc.time() - tx #244.36 

colnames(luad.clara$medoids) <- rownames(luad.2)

#Model survival outcomes
library(survival)
Surv.luad <- Surv(event=luad.surv.2$X_EVENT, time=luad.surv.2$X_TIME_TO_EVENT)

library(glmnet)
luadCOX.cv <- cv.glmnet(t(luad.clara$medoids), Surv.luad, family="cox", alpha=1, nfolds=10)
luadCOX <- glmnet(t(luad.clara$medoids), Surv.luad, family="cox", alpha=1, lambda=luadCOX.cv$lambda.min)

library(randomForestSRC)
temp <- cbind.data.frame("event"=luad.surv.2$X_EVENT, "time"=luad.surv.2$X_TIME_TO_EVENT, t(luad.clara$medoids))
luadRF <- rfsrc(Surv(time, event) ~ ., data=temp, ntree = 500, block.size = 10, importance=TRUE)

par(cex=0.25, mar=c(1,14,2,2))
barplot(sort(luadRF$importance), horiz=TRUE, las=1)

# other clinical variables
x <- match(colnames(luad.clara$medoids), luad.pheno$submitter_id.samples)
luad.pheno.2 <- luad.pheno[x,]
luad.pheno.factor <- apply(luad.pheno.2, 2, as.factor)
#

temp2 <- cbind.data.frame("event"=luad.surv.2$X_EVENT, "time"=luad.surv.2$X_TIME_TO_EVENT, luad.pheno.factor)
luadRF2 <- rfsrc(Surv(time, event) ~ ., data=temp2, ntree = 500, block.size = 10, na.action="na.impute")

par(cex=0.25, mar=c(1,14,2,2))
barplot(sort(luadRF2$importance), horiz=TRUE, las=1)

temp3 <- cbind.data.frame("event"=luad.surv.2$X_EVENT, "time"=luad.surv.2$X_TIME_TO_EVENT, luad.pheno.factor, t(luad.clara$medoids))
luadRF3 <- rfsrc(Surv(time, event) ~ ., data=temp3, ntree = 500, block.size = 10, na.action="na.impute")


temp4 <- cbind.data.frame("event"=luad.surv.2$X_EVENT, "time"=luad.surv.2$X_TIME_TO_EVENT, luad.2)
luadRF4 <- rfsrc(Surv(time, event) ~ ., data=temp4, ntree = 500, block.size = 10, na.action="na.impute")

temp5 <- cbind.data.frame("nt"=luad.pheno.2$new_tumor_event_after_initial_treatment, t(luad.clara$medoids))
luadRF5 <- rfsrc(as.factor(nt) ~ ., data=temp5, ntree = 500, block.size = 10, na.action="na.omit")

temp6 <- cbind.data.frame("nt"=luad.pheno.2$new_tumor_event_after_initial_treatment, luad.2)
luadRF6 <- rfsrc(as.factor(nt) ~ ., data=temp6, ntree = 500, block.size = 10, na.action="na.omit")


#testfit1 <- coxph(Surv(healthvalue, rep(1,nrow(y) ) ) ~numberofdrugs+treatment+improved, data= y)
#library(car) 
#cvif <- vif(testfit1) #assumes testfit from :  lrm, ols, psm, cph, Rq, Glm, glm


library(caret)
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
results <- rfe(t(luad.clara$medoids), Surv(time=luad.surv.2$X_TIME_TO_EVENT, event=luad.surv.2$X_EVENT), 
               rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))