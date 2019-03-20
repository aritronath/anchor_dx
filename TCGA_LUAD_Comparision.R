load("C:/Users/nath0050/Desktop/TCGA_lung/TCGA_lung.RData")

#### Clean-up ####
var.filter <-  function (dat, cutoff=0) {
  vars <- apply(dat, 1, var)
  dat <- dat[-which(vars == cutoff), ]
  return(dat)
}

temp <- var.filter(luad)

library(glmnet)



#testfit1 <- coxph(Surv(healthvalue, rep(1,nrow(y) ) ) ~numberofdrugs+treatment+improved, data= y)
#library(car) 
#cvif <- vif(testfit1) #assumes testfit from :  lrm, ols, psm, cph, Rq, Glm, glm

