#### Match cluster assignment between GDSC and CTRP/CCLE RNA ####
y <- match(colnames(rna), colnames(ctrp.rna))
new.gex <- ctrp.rna[, na.omit(y)]

new.clg <- clg_50[which(!is.na(y))]
pc50.ctrp <- get_pc1(new.gex, new.clg)
new.clg <- clg_52.5[which(!is.na(y))]
pc52.5.ctrp <- get_pc1(new.gex, new.clg)
new.clg <- clg_55[which(!is.na(y))]
pc55.ctrp <- get_pc1(new.gex, new.clg)

library(glmnet)
#all genes
ctrp_IKres <- list()
for (i in 1:ncol(ctrp.auc)) {
  ax <- which(!is.na(ctrp.auc[,i]))
  auc <- ctrp.auc[ax, i]
  x <- match (names(auc), rownames(ctrp.rna))
  rexp <- ctrp.rna[x,] 
  
  ctrp_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}

#Tree-cut 50 
pc50.ctrp_IKres <- list()
rownames(pc50.ctrp) <- rownames(new.gex)
for (i in 1:ncol(ctrp.auc)) {
  ax <- which(!is.na(ctrp.auc[,i]))
  auc <- ctrp.auc[ax, i]
  x <- match (names(auc), rownames(new.gex))
  rexp <- pc50.ctrp[x,] 
  
  pc50.ctrp_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}

#Tree-cut 52.5
pc52.5.ctrp_IKres <- list()
rownames(pc52.5.ctrp) <- rownames(new.gex)
for (i in 1:ncol(ctrp.auc)) {
  ax <- which(!is.na(ctrp.auc[,i]))
  auc <- ctrp.auc[ax, i]
  x <- match (names(auc), rownames(new.gex))
  rexp <- pc52.5.ctrp[x,] 
  
  pc52.5.ctrp_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}

#Tree-cut 55
pc55.ctrp_IKres <- list()
rownames(pc55.ctrp) <- rownames(new.gex)
for (i in 1:ncol(ctrp.auc)) {
  ax <- which(!is.na(ctrp.auc[,i]))
  auc <- ctrp.auc[ax, i]
  x <- match(names(auc), rownames(pc55.ctrp))
  rexp <- pc55.ctrp[x,]
  
  pc55.ctrp_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}

pc55.ctrp_IKr <- unlist(sapply(pc55.ctrp_IKres, '[', 3))
pc52.5.ctrp_IKr <- unlist(sapply(pc52.5.ctrp_IKres, '[', 3))
pc50.ctrp_IKr <- unlist(sapply(pc50.ctrp_IKres, '[', 3))
ctrp_IKr <- unlist(sapply(ctrp_IKres, '[', 3))

boxplot(pc55.ctrp_IKr, pc52.5.ctrp_IKr, pc50.ctrp_IKr, ctrp_IKr)

t.test(pc55.ctrp_IKr, ctrp_IKr)
t.test(pc52.5.ctrp_IKr, ctrp_IKr)
t.test(pc50.ctrp_IKr, ctrp_IKr)

# CTRP medoids 90
clara90.ctrp_IKres <- list()
for (i in 1:ncol(ctrp.auc)) {
  ax <- which(!is.na(ctrp.auc[,i]))
  auc <- ctrp.auc[ax, i]
  
  x <- match(names(auc), rownames(ctrp.rna))
  y <- match(ctrp.90med, colnames(ctrp.rna))
  rexp <- ctrp.rna[x,y]
  
  clara90.ctrp_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}

clara90.ctrp_IKr <- unlist(sapply(clara90.ctrp_IKres, '[', 3))
boxplot(pc55.ctrp_IKr, pc52.5.ctrp_IKr, pc50.ctrp_IKr, clara90.ctrp_IKr, ctrp_IKr)

### Random sampling for prediction ####
n1 <- ncol(pc55); n2 <- ncol(pc52.5); n3 <- ncol(pc50); 

set.seed(1001)
cx1 <- sample(ncol(ctrp.rna), n1)

set.seed(1001)
cx2 <- sample(ncol(ctrp.rna), n2)

set.seed(1001)
cx3 <- sample(ncol(ctrp.rna), n3)

cx1_IKres <- list(); cx2_IKres <- list(); cx3_IKres <- list(); 

for (i in 1:ncol(ctrp.auc)) {
  ax <- which(!is.na(ctrp.auc[,i]))
  auc <- ctrp.auc[ax,i]
  
  cx1_IKres[[i]] <- kfoldCV(ctrp.rna[ax,cx1], auc, k=5, a=0, int=TRUE, sr=FALSE)
  cx2_IKres[[i]] <- kfoldCV(ctrp.rna[ax,cx2], auc, k=5, a=0, int=TRUE, sr=FALSE)
  cx3_IKres[[i]] <- kfoldCV(ctrp.rna[ax,cx3], auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}

cx1_IKr <- unlist(sapply(cx1_IKres, '[', 3))
cx2_IKr <- unlist(sapply(cx2_IKres, '[', 3))
cx3_IKr <- unlist(sapply(cx3_IKres, '[', 3))

boxplot(cx1_IKr, pc55.ctrp_IKr, cx2_IKr, pc52.5.ctrp_IKr, cx3_IKr, pc50.ctrp_IKr, clara90.ctrp_IKr, ctrp_IKr)

save(cx1_IKres, cx2_IKres, cx3_IKres, 
     pc50.ctrp_IKres, pc52.5.ctrp_IKres, pc55.ctrp_IKres, 
     clara90.ctrp_IKres, ctrp_IKres, file="CTRP_PC_CX_IKres.RData")

####How similar are eigen genes between two datasets? ####
temp <- toupper(gsub("-", "", rownames(new.gex)))
temp <- gsub("\\.", "", temp)
x <- match(rownames(pc50), temp)

cor.test(pc50[which(!is.na(x)), 4], pc50.ctrp[na.omit(x), 5], method="spearman")

