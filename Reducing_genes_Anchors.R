### Indentify a "medoid" gene for the cluster ###
#different methods
y <- which(clg_50==1)
temp <- auc.gdsc.rna[,y]

#center - mean expression of the cluster 
mean_temp <- apply(temp, 1, mean)

#based on correlation with the mean
cor_temp <- apply(temp, 2, function (x) cor(x, mean_temp))
which(max(cor_temp)==cor_temp)
which(cor_temp==max(cor_temp))
plot(mean_temp, temp[,225])

#based on euclidean distance from the mean 
dist_temp <- apply(temp, 2, function (x) dist(t(cbind(x, mean_temp))))
which(dist_temp==min(dist_temp))
plot(mean_temp, temp[,482])

#compared to PC1 of the cluster 
x <- match(rownames(temp), rownames(pc50))
plot(pc50[x,1], temp[,225])
plot(pc50[x,1], temp[,482])

plot(pc50[x,1], mean_temp) #pc1 and mean are exactly the same!

cor(pc50[x,1], temp[,225])
cor(pc50[x,1], temp[,482])

#### Predicting sensitivity using cluster proxies #### 
#function to identify proxies from clusters = "anchors"
get_anchor <- function(clust, exp) {
  exp <- exp[!duplicated(rownames(exp)),]
  anchor_genes <- data.frame(matrix(nrow=nrow(exp), ncol=length(unique(clust))), row.names=rownames(exp))
  
  for (i in 1:length(unique(clust))) {
    y <- which(clust==i)
    temp <- exp[,y]
    mean_temp <- apply(temp, 1, mean)
    cor_temp <- apply(temp, 2, function (x) cor(x, mean_temp))
    agi <- which(cor_temp==max(abs(cor_temp)))
    anchor_genes[,i] <- temp[,agi]
    colnames(anchor_genes)[i] <- colnames(temp)[agi]
  }
  return(anchor_genes)
}

pc50_anchors <- get_anchor(clg_50, auc.gdsc.rna)
pc52.5_anchors <- get_anchor(clg_52.5, auc.gdsc.rna)
pc55_anchors <- get_anchor(clg_55, auc.gdsc.rna)

get_pc_anchor <- function (clust, pcmat, exp) {
  anchor_genes <- matrix(data=NA, nrow=nrow(exp), ncol=ncol(pcmat))
  rownames(anchor_genes) <- rownames(exp)
  
  for (i in 1:ncol(pcmat)) {
    y <- which(clust==i)
    temp <- exp[,y]
    cor_temp <- apply(temp, 2, function (x) cor(x, pcmat[,i]))
    agi <- which(abs(cor_temp)==max(abs(cor_temp)))
    anchor_genes[,i] <- temp[,agi]
  }
  return(anchor_genes)
}

pc50_pc_anchors <- get_pc_anchor(clg_50, pc50, rna)


### Sensitivity prediction ###
pc55_anchors_IKres <- list()
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc55_anchors))
  rexp <- pc55_anchors[x,]
  pc55_anchors_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  print (i)  
}

pc52.5_anchors_IKres <- list()
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc52.5_anchors))
  rexp <- pc52.5_anchors[x,]
  pc52.5_anchors_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  print (i)  
}

pc50_anchors_IKres <- list()
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc50_anchors))
  rexp <- pc50_anchors[x,]
  pc50_anchors_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  print (i)  
}

pc50_pc_anchors_IKres <- list()
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(pc50_pc_anchors))
  rexp <- pc50_pc_anchors[x,]
  pc50_pc_anchors_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  print (i)  
}


pc55_anchors_IKr <- unlist(sapply(pc55_anchors_IKres, '[', 3))
pc52.5_anchors_IKr <- unlist(sapply(pc52.5_anchors_IKres, '[', 3))
pc50_anchors_IKr <- unlist(sapply(pc50_anchors_IKres, '[', 3))
pc50_pc_anchors_IKr <- unlist(sapply(pc50_pc_anchors_IKres, '[', 3))

boxplot(pc55_IKr, pc55_anchors_IKr, pc52.5_IKr, pc52.5_anchors_IKr, pc50_IKr, pc50_anchors_IKr, pc50_pc_anchors_IKr)
boxplot(sx1_IKr, pc55_anchors_IKr, sx2_IKr, pc52.5_anchors_IKr, sx3_IKr, pc50_anchors_IKr, pc50_pc_anchors_IKr)

y <- match(rownames(auc.gdsc.rna), rownames(pc50))
plot(pc50[y,9], pc50_anchors[,9])

#### Predicting sensitivity using PC proxies ####
get_IPC_anchor <- function (pcmat, exp) {
  anchor_genes <- matrix(data=NA, nrow=nrow(exp), ncol=ncol(pcmat))
  rownames(anchor_genes) <- rownames(exp)
  
  for (i in 1:ncol(pcmat)) {
    cor_temp <- apply(exp, 2, function (x) cor(x, pcmat[,i]))
    agi <- which(abs(cor_temp)==max(abs(cor_temp)))
    anchor_genes[,i] <- exp[,agi]
    message(print(i))
  }
  return(anchor_genes)
}

IPC15_anchors <- get_IPC_anchor(gdsc.pc15, auc.gdsc.rna)
IPC10_anchors <- get_IPC_anchor(gdsc.pc10, auc.gdsc.rna)
IPC5_anchors <- get_IPC_anchor(gdsc.pc5, auc.gdsc.rna)

IPC15_anchors_IKres <- list()
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(IPC15_anchors))
  rexp <- IPC15_anchors[x,]
  IPC15_anchors_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  print (i)  
}

IPC10_anchors_IKres <- list()
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(IPC10_anchors))
  rexp <- IPC10_anchors[x,]
  IPC10_anchors_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  print (i)  
}

IPC5_anchors_IKres <- list()
for (i in 1:ncol(auc.gdsc)) {
  auc <- auc.gdsc[which(!is.na(auc.gdsc[,i])), i]
  x <- match(names(auc), rownames(IPC5_anchors))
  rexp <- IPC5_anchors[x,]
  IPC5_anchors_IKres[[i]] <- kfoldCV(rexp, auc, k=5, a=0, int=TRUE, sr=FALSE)
  print (i)  
}


IPC15_anchors_IKr <- unlist(sapply(IPC15_anchors_IKres, '[', 3))
IPC10_anchors_IKr <- unlist(sapply(IPC10_anchors_IKres, '[', 3))
IPC5_anchors_IKr <- unlist(sapply(IPC5_anchors_IKres, '[', 3))

boxplot(IPC5_anchors_IKr, pc55_anchors_IKr, IPC10_anchors_IKr, pc52.5_anchors_IKr, IPC15_anchors_IKr, pc50_anchors_IKr)
boxplot(sx3_IKr, pc50_anchors_IKr, pc50_pc_anchors_IKr, IPC10_anchors_IKr)
