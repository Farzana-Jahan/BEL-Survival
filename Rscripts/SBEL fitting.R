## Fitting model to survival data using BELSpatial package and SBEL models
library(tidyverse)
library(plyr)
#library(spdep)
library(gmm)
library(emplik)
library(BELSpatial)
ospg.model.data<- read_csv("Data/ospg.data.csv")
W<- readRDS("Data/W.RDS")
ni<-rowSums(W) # no. of neighbours for each area
R<-diag(ni)
for(i in 1:nrow(R))
{
  R[i,which(W[i,]==1)]<- -1
}

ospg.model.data<-mutate(ospg.model.data,cities=ifelse(region==1,1,0),
                        regional=ifelse(region==2,1,0),
                        remote=ifelse(region==3,1,0))
y<- ospg.model.data$y.ospg.male
n<-length(y)
x<-cbind(ospg.model.data$cities,ospg.model.data$regional,ospg.model.data$remote)
p<- dim(x)[2] # no. of covariates
alpha_1<-1 # hyperparamter for tau prior
alpha_2<-0.01 # hyperparamter for tau prior
tau_inv_init<- rgamma(1,alpha_1,alpha_2) # using IG prior(1,1) for tau_inv
tau_init<- 1/tau_inv_init
g<- 10# G prior evaluated at 10 for regression coefficients' prior (Zellner prior)
prior_mean_beta<- rep(0,p) # p is the number of regression parameters, in case of one covariate, p=2
beta_init<- rnorm(3,prior_mean_beta, (1/g)*tau_inv_init)
wi_init<- 1/length(y) # y be the response variable from the data
psi_init <- rep(0,n)
var<- ospg.model.data$sd.ospg.male
# calculating MELE of Beta, beta_mele
wi=wi_init
# using gmm package to calculate initial values of beta
g<- y~ x[,1]+x[,2]+x[,3]-1
H<-x
beta_mele<- unname(gel(g,H,c(0,0,0))$coefficients)
mu_init<- exp(x%*% beta_mele + psi_init)
beta_init<-beta_mele
wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi_mu<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 
wi<-wi_mu
# SBEL ind
# fitting BEL BYM model taking rho= 1
library(parallel)
cluster<-makeCluster(3)
#clusterEvalQ(cl=cluster,.libPaths("c:/software/Rpackages"))
clusterEvalQ(cl=cluster,library(BELSpatial))
clusterExport(cl=cluster,varlist = c("y","x","n","p","var","beta_init", "psi_init", "tau_init","R", "wi"))
SBEL_BYM_ospg_male_surv<-clusterApply(cl=cluster, x=1:3, function(z){BEL_leroux_new(y,x,n,p,var,rho=0,niter=1000,
                                                                                 beta_init, psi_init, tau_init,R, wi, sd_psi=0.006, 
                                                                                 sd_beta=0.0009, sd_tau=0.4)})
save(SBEL_BYM_ospg_male_surv,file="Results/SBEL_BYM_ospg_male_surv_1000.RData")



