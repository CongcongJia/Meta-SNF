

library(MASS)
set.seed(123456)
n=200 #number of patients
p=1000 #number of features
info_num=100 #number of signal features used for clustering
noise_num=p-info_num #number of noise features
ns=50 # number of patients in each cluster
C=4 #number of clusters
truelabel=c(rep(1,ns),rep(2,ns),rep(3,ns),rep(4,ns))
value_var=6

##################################################
#Generate data   
data1=matrix(0,n,p)
data1[51:150,1:info_num]=1
data1[151:200,1:info_num]=3
data1=data1+mvrnorm(n,mu=rep(0,p),Sigma =diag(value_var,nrow=p,ncol=p))

data2=matrix(0,n,p)
data2[51:100,1:info_num]=2
data2[151:200,1:info_num]=3
data2=data2+mvrnorm(n,mu=rep(0,p),Sigma =diag(value_var,nrow=p,ncol=p))

data3=matrix(0,n,p)
data3[1:100,1:info_num]=2
data3[101:150,1:info_num]=1
data3[151:200,1:info_num]=3
data3=data3+mvrnorm(n,mu=rep(0,p),Sigma =diag(value_var,nrow=p,ncol=p))


