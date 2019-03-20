#### PCR/iPC method for prediction ####
x <- match(rownames(auc.gdsc), rownames(rna))
auc.gdsc.rna <- rna[x,]  
gdsc.pc15 <- svd(auc.gdsc.rna, nu=15, nv=0)$u
gdsc.pc10 <- svd(auc.gdsc.rna, nu=10, nv=0)$u
gdsc.pc5 <- svd(auc.gdsc.rna, nu=5, nv=0)$u

iPC15_IKres <- list(); iPC10_IKres <- list(); iPC5_IKres <- list()
rownames(gdsc.pc15) <- rownames(auc.gdsc.rna)
rownames(gdsc.pc10) <- rownames(auc.gdsc.rna)
rownames(gdsc.pc5) <- rownames(auc.gdsc.rna)

for (i in 1:ncol(auc.gdsc)) {
  ax <- which(!is.na(auc.gdsc[,i]))
  auc <- auc.gdsc[ax, i]
  
  iPC15_IKres[[i]] <- kfoldCV(gdsc.pc15[ax,], auc, k=5, a=0, int=TRUE, sr=FALSE)
  iPC10_IKres[[i]] <- kfoldCV(gdsc.pc10[ax,], auc, k=5, a=0, int=TRUE, sr=FALSE)
  iPC5_IKres[[i]] <- kfoldCV(gdsc.pc5[ax,], auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}

save(iPC15_IKres, iPC10_IKres, iPC5_IKres, file="iPC_IKres.RData")

iPC15_IKr <- unlist(sapply(iPC15_IKres, '[', 3))
iPC10_IKr <- unlist(sapply(iPC10_IKres, '[', 3))
iPC5_IKr <- unlist(sapply(iPC5_IKres, '[', 3))

boxplot(iPC15_IKr, iPC10_IKr, iPC5_IKr, pc55_IKr, pc42.5_IKr)

#### Random sampling for prediction ####
n1 <- ncol(pc55); n2 <- ncol(pc52.5)
n3 <- ncol(pc50); n4 <- ncol(pc47.5)
n5 <- ncol(pc45); n6 <- ncol(pc42.5)

set.seed(1001)
sx1 <- sample(ncol(auc.gdsc.rna), n1)

set.seed(1001)
sx2 <- sample(ncol(auc.gdsc.rna), n2)

set.seed(1001)
sx3 <- sample(ncol(auc.gdsc.rna), n3)

set.seed(1001)
sx4 <- sample(ncol(auc.gdsc.rna), n4)

set.seed(1001)
sx5 <- sample(ncol(auc.gdsc.rna), n5)

set.seed(1001)
sx6 <- sample(ncol(auc.gdsc.rna), n6)

sx1_IKres <- list(); sx2_IKres <- list(); sx3_IKres <- list(); 
sx4_IKres <- list(); sx5_IKres <- list(); sx6_IKres <- list(); 

sx1.exp <- auc.gdsc.rna[,sx1]; sx2.exp <- auc.gdsc.rna[,sx2]; sx3.exp <- auc.gdsc.rna[,sx3]
sx4.exp <- auc.gdsc.rna[,sx4]; sx5.exp <- auc.gdsc.rna[,sx5]; sx6.exp <- auc.gdsc.rna[,sx6]

for (i in 1:ncol(auc.gdsc)) {
  ax <- which(!is.na(auc.gdsc[,i]))
  auc <- auc.gdsc[ax, i]
  
  sx1_IKres[[i]] <- kfoldCV(sx1.exp[ax,], auc, k=5, a=0, int=TRUE, sr=FALSE)
  sx2_IKres[[i]] <- kfoldCV(sx2.exp[ax,], auc, k=5, a=0, int=TRUE, sr=FALSE)
  sx3_IKres[[i]] <- kfoldCV(sx3.exp[ax,], auc, k=5, a=0, int=TRUE, sr=FALSE)
  sx4_IKres[[i]] <- kfoldCV(sx4.exp[ax,], auc, k=5, a=0, int=TRUE, sr=FALSE)
  sx5_IKres[[i]] <- kfoldCV(sx5.exp[ax,], auc, k=5, a=0, int=TRUE, sr=FALSE)
  sx6_IKres[[i]] <- kfoldCV(sx6.exp[ax,], auc, k=5, a=0, int=TRUE, sr=FALSE)
  
  print (i)  
}

sx1_IKr <- unlist(sapply(sx1_IKres, '[', 3)); sx2_IKr <- unlist(sapply(sx2_IKres, '[', 3))
sx3_IKr <- unlist(sapply(sx3_IKres, '[', 3)); sx4_IKr <- unlist(sapply(sx4_IKres, '[', 3))
sx5_IKr <- unlist(sapply(sx5_IKres, '[', 3)); sx6_IKr <- unlist(sapply(sx6_IKres, '[', 3))

boxplot(sx1_IKr, pc55_IKr, sx2_IKr,  pc52.5_IKr, sx3_IKr, pc50_IKr, sx4_IKr, pc47.5_IKr, sx5_IKr, pc45_IKr, 
        sx6_IKr, pc42.5_IKr, notch=TRUE, col=c("red", "gold"))


  #### FIGURE: Compare prediction accuracy cluster vs. random sampled genes ####
library(reshape2)
library(ggplot2)
temp <- melt(cbind(sx1_IKr, pc55_IKr, sx2_IKr,  pc52.5_IKr, sx3_IKr, pc50_IKr, sx4_IKr, pc47.5_IKr, 
                   sx5_IKr, pc45_IKr, sx6_IKr, pc42.5_IKr))
grp <- rep(c(rep("Random", 265), rep("Cluster", 265)), 6)
ssize = c(rep("N=7", 530), rep("N=22", 530), rep("N=90", 530), rep("N=268", 530), 
          rep("N=863", 530), rep("N=2106", 530))
temp <- cbind(temp, grp, ssize)

a <- ggplot(temp, aes(x=reorder(ssize, value, FUN=median), y=value, fill=reorder(grp, value, FUN=median, .desc=FALSE))) 
a + geom_boxplot(notch=TRUE, width=0.5, alpha=0.8) +
  labs(x="Number of training features", y="Average Spearman correlation",  fill="Model") + 
  scale_fill_manual(values=c("skyblue", "red")) + theme_classic(base_size = 16)
  
  #Annotation 
t.test(sx1_IKr, pc55_IKr, paired=T) #P < 2.2e-16
t.test(sx2_IKr, pc52.5_IKr, paired=T) #P < 2.2e-16
t.test(sx3_IKr, pc50_IKr, paired=T) #P = 1.61e-8
t.test(sx4_IKr, pc47.5_IKr, paired=T) #P = 0.00058
t.test(sx5_IKr, pc45_IKr, paired=T) #P = 0.376
t.test(sx6_IKr, pc42.5_IKr, paired=T) #P = 0.117

  #####

  ##### Figure: comparing order of drugs in each analysis ####
par(mar=c(3,4.5,1,2), mfrow=c(2,3), cex=1.1, lwd=2, pch=16, bty="l")
plot(pc55_IKr~gdsc_IKr, ylab="Cluster N=7", xlab=""); abline(lm(pc55_IKr~gdsc_IKr))
plot(pc52.5_IKr~gdsc_IKr, ylab="Cluster N=22", xlab=""); abline(lm(pc52.5_IKr~gdsc_IKr))
plot(pc50_IKr~gdsc_IKr, ylab="Cluster N=90", xlab=""); abline(lm(pc50_IKr~gdsc_IKr))
plot(pc47.5_IKr~gdsc_IKr, ylab="Cluster N=268", xlab=""); abline(lm(pc47.5_IKr~gdsc_IKr))
plot(pc45_IKr~gdsc_IKr, ylab="Cluster N=863", xlab=""); abline(lm(pc45_IKr~gdsc_IKr))
plot(pc42.5_IKr~gdsc_IKr, ylab="Cluster N=2106", xlab=""); abline(lm(pc42.5_IKr~gdsc_IKr))

cor(pc55_IKr,gdsc_IKr, method="spearman") #rho=0.83
cor(pc52.5_IKr,gdsc_IKr, method="spearman") #rho=0.89
cor(pc50_IKr,gdsc_IKr, method="spearman") #rho=0.94
cor(pc47.5_IKr,gdsc_IKr, method="spearman") #rho=0.96
cor(pc45_IKr,gdsc_IKr, method="spearman") #rho=0.97
cor(pc42.5_IKr,gdsc_IKr, method="spearman") #rho=0.98

  #####

#### How many variables are needed for well-predicted drugs ####
intersect(which(gdsc_IKr > 0.5), which(pc55_IKr > 0.5))

#### What do these predictors predict? ####
### CFEs? Tissue type? Drug target? ###
auc.CFE
auc.tiss






