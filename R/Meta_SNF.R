
###############load Packages
library(NMF)
library(SNFtool)
library(survival)
library(survminer)
library(survRM2)
library(ggplot2)

###############define a seed and set all the parameters
gg <- 123456
K = 20      ##number of neighbors, usually (10~30)
alpha = 0.5 ##hyperparameter, usually (0.3~0.8)
T = 20      ###Number of iterations, usually (10~50)


###############load real datasets
setwd("~/Meta-SNF/KIRP/data")
load("methy.Rdata")
load("miRNA.Rdata")
load("mRNA.Rdata")
load("clinic_KIRP_match.Rdata")

survdat=cbind(clinic_KIRP_match[,1],as.numeric(clinic_KIRP_match[,7]),as.numeric(clinic_KIRP_match[,2]))
rownames(survdat)=survdat[,1]
survdat=survdat[,-1]
colnames(survdat)=c("time","event")
survdat=apply(survdat,2,as.numeric)
survdat=as.data.frame(survdat)

######################
# run Meta-SNF
######################

###############Step1 Extract the metagene matrices by NMF
###############estimate rank for datasets
rank.best<-function(coph){
  coph_diff <- NULL
  for (i in 2:length(coph)) 
  {
    coph_diff <- c(coph_diff, coph[i-1]-coph[i])
  }
  k.best <- which.max(coph_diff)+1
  return(k.best)
}

est_miRNA<- nmf(miRNA, 2:11, nrun = 30,seed = gg)
coph_miRNA<-est_miRNA[["measures"]][["cophenetic"]]
k.best_miRNA<-rank.best(coph_miRNA)
est_mRNA <- nmf(mRNA, 2:11, nrun = 30,seed = gg)
coph_mRNA<-est_mRNA[["measures"]][["cophenetic"]]
k.best_mRNA<-rank.best(coph_mRNA)
est_methy <- nmf(methy, 2:11, nrun = 30,seed = gg)
coph_methy<-est_methy[["measures"]][["cophenetic"]]
k.best_methy<-rank.best(coph_methy)


result_mRNA<-nmf(mRNA,k.best_mRNA,nrun=30)
result_miRNA<-nmf(miRNA,k.best_miRNA,nrun=30)
result_methy<-nmf(methy,k.best_methy,nrun=30)

H_mRNA<- result_mRNA@fit@H
H_miRNA <- result_miRNA@fit@H
H_methy <- result_methy@fit@H

###############Step2 Integrate metagene matrices by SNF
Dist1=(dist2(as.matrix(t(H_mRNA)),as.matrix(t(H_mRNA))))^(1/2)
Dist2=(dist2(as.matrix(t(H_miRNA)),as.matrix(t(H_miRNA))))^(1/2)
Dist3=(dist2(as.matrix(t(H_methy)),as.matrix(t(H_methy))))^(1/2)

W1<-affinityMatrix(Dist1,K, alpha)
W2<-affinityMatrix(Dist2,K, alpha)
W3<-affinityMatrix(Dist3,K, alpha)
W_Meta= SNF(list(W1,W2,W3), K, T)

###############estimates the number of clusters
estimationResult_NMF = estimateNumberOfClustersGivenGraph(W_Meta, NUMC=3:5);
num_NMF_gap=estimationResult_NMF$`Eigen-gap best` 
num_NMF_cost=estimationResult_NMF$`Rotation cost best`

###############apply the spectral clustering method
group = spectralClustering(W_Meta,num_NMF_gap)
table(group)

###############Plot survival curves
surv=cbind(survdat,group)
sdf=survdiff(Surv(survdat[,1],as.numeric(survdat[,2]))~surv$group,data=surv)
p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)

fit=survfit(Surv(survdat[,1],as.numeric(survdat[,2]))~surv$group,data=surv)
windowsFonts(myFont = windowsFont("Times New Roman")) 
tiff(filename = "Meta-SNF.tiff",res = 300,width = 4.5,height = 2.5,units = 'in',compression = c("lzw"))
ggsurvplot(fit,conf.int=F,xlab='Survival time(months)',
           risk.table =F, 
           palette = c("#377EB8", "#E41A1C","#4DAF4A","#FDAE61"),
           legend.title = "Subtype", 
           legend.labs = c("Cluster1", "Cluster2","Cluster3","Cluster4"),
           legend = "right",
           font.x = 14 ,font.y = 14 ,font.tickslab=12,font.legend=14,
           linetype = 1,size=1.3)
dev.off()