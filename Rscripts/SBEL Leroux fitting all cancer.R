## Fitting model to survival data using BELSpatial package and SBEL models
library(tidyverse)
library(plyr)
#library(spdep)
library(gmm)
library(emplik)
library(BELSpatial)
all.model.data<- read_csv("Data/all.data.csv")
cov<- read_csv("Data/SA2 Data_v2.csv")
cov<-mutate(cov,region=ra5cat) 
cov$region[cov$region==3]<-2
cov$region[cov$region==4]<-3
cov$region[cov$region==5]<-3
cov$region<- factor(cov$region)
summary(cov$region)
head(cov)
library(dplyr)
cov<-rename(cov,c("sa2"="SA2_code"))
all.model.data<- inner_join(all.model.data,cov,by="SA2_code")
W<- readRDS("Data/W.RDS")
ni<-rowSums(W) # no. of neighbours for each area
R<-diag(ni)
for(i in 1:nrow(R))
{
  R[i,which(W[i,]==1)]<- -1
}

all.model.data<-mutate(all.model.data,cities=ifelse(region==1,1,0),
                        regional=ifelse(region==2,1,0),
                        remote=ifelse(region==3,1,0))
y<- all.model.data$y.all.male
n<-length(y)
x<-cbind(all.model.data$cities,all.model.data$regional,all.model.data$remote)
p<- dim(x)[2] # no. of covariates

# changing initial value from last fit
alpha_1<-1 # hyperparamter for tau prior
alpha_2<-0.01 # hyperparamter for tau prior
tau_inv_init<- rgamma(1,alpha_1,alpha_2) # using IG prior(1,1) for tau_inv
tau_init<- 1/tau_inv_init
g<- 10# G prior evaluated at 10 for regression coefficients' prior (Zellner prior)
prior_mean_beta<- rep(0,p) # p is the number of regression parameters, in case of one covariate, p=2
beta_init<- rnorm(3,prior_mean_beta, (1/g)*tau_inv_init)
wi_init<- 1/length(y) # y be the response variable from the data
psi_init <- rep(0,n)
var<- all.model.data$sd.all.male^2
# calculating MELE of Beta, beta_mele
wi=wi_init
# using gmm package to calculate initial values of beta
g<- y~ x[,1]+x[,2]+x[,3]-1
H<-x
beta_mele<- unname(gel(g,H,c(0,0,0))$coefficients)
mu_init<- x%*% beta_mele + psi_init
beta_init<-beta_mele
wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi_mu<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 
wi<-wi_mu
# SBEL Leroux
# fitting BEL Leroux model taking 
library(parallel)
cluster<-makeCluster(3)
#clusterEvalQ(cl=cluster,.libPaths("c:/software/Rpackages"))
clusterEvalQ(cl=cluster,library(BELSpatial))
clusterExport(cl=cluster,varlist = c("y","x","n","p","var","beta_init", "psi_init", "tau_init","R", "wi"))
SBEL_Leroux_all_male_surv<-clusterApply(cl=cluster, x=1:3, function(z){BEL_leroux_new(y,x,n,p,var,rho=0.9,niter=1000,
                                                                                      beta_init, psi_init, tau_init,R, wi, sd_psi=0.001, 
                                                                                      sd_beta=0.0005, sd_tau=1.3)})
#Checkpoint 1
save(SBEL_Leroux_all_male_surv,file="Results/SBEL_Leroux_all_male_surv_1000.RData")

# changing initial value from last fit
# using IG prior(1,1) for tau_inv
tau_init<- mean(c(SBEL_Leroux_all_male_surv[[1]]$tau[1000],SBEL_Leroux_all_male_surv[[2]]$tau[1000],
                  SBEL_Leroux_all_male_surv[[3]]$tau[1000]))
g<- 10# G prior evaluated at 10 for regression coefficients' prior (Zellner prior)
prior_mean_beta<- rep(0,p) # p is the number of regression parameters, in case of one covariate, p=2
beta_init<- colMeans(matrix(c(SBEL_Leroux_all_male_surv[[1]]$Beta[,1000],SBEL_Leroux_all_male_surv[[2]]$Beta[,1000],
                              SBEL_Leroux_all_male_surv[[3]]$Beta[,1000]),nrow=3,byrow = F))
# y be the response variable from the data
psi_init <- colMeans(matrix(c(SBEL_Leroux_all_male_surv[[1]]$psi[,1000],SBEL_Leroux_all_male_surv[[2]]$psi[,1000],
                              SBEL_Leroux_all_male_surv[[3]]$psi[,1000]),nrow=2148,byrow = F))
var<- all.model.data$sd.all.male^2
# calculating MELE of Beta, beta_mele

# using gmm package to calculate initial values of beta
mu_init<- x%*% beta_init+ psi_init
wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 

# SBEL Leroux
# fitting BEL Leroux model taking rho= 1
library(parallel)
cluster<-makeCluster(3)
#clusterEvalQ(cl=cluster,.libPaths("c:/software/Rpackages"))
clusterEvalQ(cl=cluster,library(BELSpatial))
clusterExport(cl=cluster,varlist = c("y","x","n","p","var","beta_init", "psi_init", "tau_init","R", "wi"))
SBEL_Leroux_all_male_surv2<-clusterApply(cl=cluster, x=1:3, function(z){BEL_leroux_new(y,x,n,p,var,rho=0.9,niter=1000,
                                                                                       beta_init, psi_init, tau_init,R, wi, sd_psi=0.001, 
                                                                                       sd_beta=0.0005, sd_tau=1.3)})
#Checkpoint 2
save(SBEL_Leroux_all_male_surv2,file="Results/SBEL_Leroux_all_male_surv_2000.RData")

# changing initial value from last fit
# using IG prior(1,1) for tau_inv
tau_init<- mean(c(SBEL_Leroux_all_male_surv2[[1]]$tau[1000],SBEL_Leroux_all_male_surv2[[2]]$tau[1000],
                  SBEL_Leroux_all_male_surv2[[3]]$tau[1000]))
g<- 10# G prior evaluated at 10 for regression coefficients' prior (Zellner prior)
prior_mean_beta<- rep(0,p) # p is the number of regression parameters, in case of one covariate, p=2
beta_init<- colMeans(matrix(c(SBEL_Leroux_all_male_surv2[[1]]$Beta[,1000],SBEL_Leroux_all_male_surv2[[2]]$Beta[,1000],
                              SBEL_Leroux_all_male_surv2[[3]]$Beta[,1000]),nrow=3,byrow = F))
# y be the response variable from the data
psi_init <- colMeans(matrix(c(SBEL_Leroux_all_male_surv2[[1]]$psi[,1000],SBEL_Leroux_all_male_surv2[[2]]$psi[,1000],
                              SBEL_Leroux_all_male_surv2[[3]]$psi[,1000]),nrow=2148,byrow = F))
var<- all.model.data$sd.all.male^2
# calculating MELE of Beta, beta_mele

# using gmm package to calculate initial values of beta
mu_init<- x%*% beta_init+ psi_init
wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 

# SBEL Leroux
# fitting BEL Leroux model taking rho= 1
library(parallel)
cluster<-makeCluster(3)
#clusterEvalQ(cl=cluster,.libPaths("c:/software/Rpackages"))
clusterEvalQ(cl=cluster,library(BELSpatial))
clusterExport(cl=cluster,varlist = c("y","x","n","p","var","beta_init", "psi_init", "tau_init","R", "wi"))
SBEL_Leroux_all_male_surv3<-clusterApply(cl=cluster, x=1:3, function(z){BEL_leroux_new(y,x,n,p,var,rho=0.9,niter=1000,
                                                                                       beta_init, psi_init, tau_init,R, wi, sd_psi=0.001, 
                                                                                       sd_beta=0.0005, sd_tau=1.3)})
#Checkpoint 3
save(SBEL_Leroux_all_male_surv3,file="Results/SBEL_Leroux_all_male_surv_3000.RData")

# changing initial value from last fit
# using IG prior(1,1) for tau_inv
tau_init<- mean(c(SBEL_Leroux_all_male_surv3[[1]]$tau[1000],SBEL_Leroux_all_male_surv3[[2]]$tau[1000],
                  SBEL_Leroux_all_male_surv3[[3]]$tau[1000]))
g<- 10# G prior evaluated at 10 for regression coefficients' prior (Zellner prior)
prior_mean_beta<- rep(0,p) # p is the number of regression parameters, in case of one covariate, p=2
beta_init<- colMeans(matrix(c(SBEL_Leroux_all_male_surv3[[1]]$Beta[,1000],SBEL_Leroux_all_male_surv3[[2]]$Beta[,1000],
                              SBEL_Leroux_all_male_surv3[[3]]$Beta[,1000]),nrow=3,byrow = F))
# y be the response variable from the data
psi_init <- colMeans(matrix(c(SBEL_Leroux_all_male_surv3[[1]]$psi[,1000],SBEL_Leroux_all_male_surv3[[2]]$psi[,1000],
                              SBEL_Leroux_all_male_surv3[[3]]$psi[,1000]),nrow=2148,byrow = F))
var<- all.model.data$sd.all.male^2
# calculating MELE of Beta, beta_mele

# using gmm package to calculate initial values of beta
mu_init<- x%*% beta_init+ psi_init
wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 

# SBEL Leroux
# fitting BEL Leroux model taking rho= 1
library(parallel)
cluster<-makeCluster(3)
#clusterEvalQ(cl=cluster,.libPaths("c:/software/Rpackages"))
clusterEvalQ(cl=cluster,library(BELSpatial))
clusterExport(cl=cluster,varlist = c("y","x","n","p","var","beta_init", "psi_init", "tau_init","R", "wi"))
SBEL_Leroux_all_male_surv4<-clusterApply(cl=cluster, x=1:3, function(z){BEL_leroux_new(y,x,n,p,var,rho=0.9,niter=1000,
                                                                                       beta_init, psi_init, tau_init,R, wi, sd_psi=0.005, 
                                                                                       sd_beta=0.001, sd_tau=1.2)})
#Checkpoint 4
save(SBEL_Leroux_all_male_surv4,file="Results/SBEL_Leroux_all_male_surv_4000.RData")

# changing initial value from last fit
# using IG prior(1,1) for tau_inv
tau_init<- mean(c(SBEL_Leroux_all_male_surv4[[1]]$tau[1000],SBEL_Leroux_all_male_surv4[[2]]$tau[1000],
                  SBEL_Leroux_all_male_surv4[[3]]$tau[1000]))
g<- 10# G prior evaluated at 10 for regression coefficients' prior (Zellner prior)
prior_mean_beta<- rep(0,p) # p is the number of regression parameters, in case of one covariate, p=2
beta_init<- colMeans(matrix(c(SBEL_Leroux_all_male_surv4[[1]]$Beta[,1000],SBEL_Leroux_all_male_surv4[[2]]$Beta[,1000],
                              SBEL_Leroux_all_male_surv4[[3]]$Beta[,1000]),nrow=3,byrow = F))
# y be the response variable from the data
psi_init <- colMeans(matrix(c(SBEL_Leroux_all_male_surv4[[1]]$psi[,1000],SBEL_Leroux_all_male_surv4[[2]]$psi[,1000],
                              SBEL_Leroux_all_male_surv4[[3]]$psi[,1000]),nrow=2148,byrow = F))
var<- all.model.data$sd.all.male^2
# calculating MELE of Beta, beta_mele

# using gmm package to calculate initial values of beta
mu_init<- x%*% beta_init+ psi_init
wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 

# SBEL Leroux
# fitting BEL Leroux model taking rho= 1
library(parallel)
cluster<-makeCluster(3)
#clusterEvalQ(cl=cluster,.libPaths("c:/software/Rpackages"))
clusterEvalQ(cl=cluster,library(BELSpatial))
clusterExport(cl=cluster,varlist = c("y","x","n","p","var","beta_init", "psi_init", "tau_init","R", "wi"))
SBEL_Leroux_all_male_surv5<-clusterApply(cl=cluster, x=1:3, function(z){BEL_leroux_new(y,x,n,p,var,rho=0.9,niter=1000,
                                                                                       beta_init, psi_init, tau_init,R, wi, sd_psi=0.001, 
                                                                                       sd_beta=0.0005, sd_tau=1.3)})
#Checkpoint 5
save(SBEL_Leroux_all_male_surv5,file="Results/SBEL_Leroux_all_male_surv_5000.RData")



