.libPaths("c:/software/Rpackages")
library(tidyverse)
library(plyr)
#library(spdep)
library(gmm)
library(CARBayes)
# exploring data
#all survival male
all.model.data<- read_csv("Data/all.data.csv")
head(all.model.data)
library(tidyverse)
ggplot(data=all.model.data, aes(x= y.all.male))+geom_density()+theme_bw()

# adding neighbourhood matrix
W<- readRDS("Data/W.RDS")
# adding covariate to the data
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
#write_csv(all.model.data,"Data/all.data.csv")
# Parametric Model fit

library(CARBayes)
all.model.data$region<-factor(all.model.data$region)
param_BYM<-S.CARleroux(y.all.male~region-1, data=all.model.data,W= W, family="gaussian",rho=1, n.sample = 100000,
                       burnin = 10000,thin=1)
save(param_BYM,file="Results/all_BYM.RData")
param_Ind<-S.CARleroux(y.all.male~region-1, data=all.model.data,W= W, family="gaussian",rho=0, n.sample = 100000,
                       burnin = 10000,thin=1)
save(param_Ind,file= "Results/all_Ind.RData")
param_Leroux<-S.CARleroux(y.all.male~region-1, data=all.model.data,W= W, family="gaussian", n.sample = 100000,
                       burnin = 10000,thin=1)
save(param_Leroux,file= "Results/all_Leroux.RData")

library(ngspatial)
param_moranBasis<-sparse.sglmm(y.all.male~region-1,data=all.model.data,A= W,
                               family="gaussian",minit = 10000, maxit=100000)
save(param_moranBasis,file= "Results/all_moranBasis.RData")

# calcuation of WAIC and formatting paramteric output
thin_param_MCMC<- function(output,burnin=10000,thin=10)
{
  if(class(output)=="numeric")
  {
    niter=length(output)
    out_postburn<-output[(burnin+1):niter]
    n<-length(out_postburn)
    out_thin<- out_postburn[seq(1,n,by=thin)]
    return(out_thin)
  }
  
  else{
    niter=dim(output)[1]
    out_postburn<-output[(burnin+1):niter,]
    n<-dim(out_postburn)[1]
    out_thin<- out_postburn[seq(1,n,by=thin),]
    return(out_thin)
  }
  out_thin
}
chain1_beta0<-as.numeric(param_Ind$samples$beta[,1])
chain1_beta1<-as.numeric(param_Ind$samples$beta[,2])
chain1_beta1<-as.numeric(param_Ind$samples$beta[,3])
chain1_psi<-param_Ind$samples$phi
chain1_tau<-as.numeric(param_Ind$samples$tau)
chain2_beta0<-as.numeric(param_Ind[[2]]$samples$beta[,1])
chain2_beta1<-as.numeric(param_Ind[[2]]$samples$beta[,2])
chain2_psi<-param_Ind[[2]]$samples$phi
chain2_tau<-as.numeric(param_Ind[[2]]$samples$tau)
chain3_beta0<-as.numeric(param_Ind[[3]]$samples$beta[,1])
chain3_beta1<-as.numeric(param_Ind[[3]]$samples$beta[,2])
chain3_psi<-param_Ind[[3]]$samples$phi
chain3_tau<-as.numeric(param_Ind[[3]]$samples$tau)

# Apply these to all BEL models after thinning and save with proper names
# thinned chains
chain1_beta0_thin<-thin_param_MCMC(chain1_beta0,burnin =0, thin = 10)
chain1_beta1_thin<-thin_param_MCMC(chain1_beta1,burnin = 0, thin = 10)
chain1_beta2_thin<-thin_param_MCMC(chain1_beta2,burnin = 0, thin = 10)
chain1_psi_thin<-thin_param_MCMC(chain1_psi,burnin = 0, thin = 10)
chain1_tau_thin<-thin_param_MCMC(chain1_tau,burnin = 0, thin = 10)
chain2_beta0_thin<-thin_param_MCMC(chain2_beta0,burnin = 0, thin = 10)
chain2_beta1_thin<-thin_param_MCMC(chain2_beta1,burnin = 0, thin = 10)
chain2_psi_thin<-thin_param_MCMC(chain2_psi,burnin = 0, thin = 10)
chain2_tau_thin<-thin_param_MCMC(chain2_tau,burnin = 0, thin = 10)
chain3_beta0_thin<-thin_param_MCMC(chain3_beta0,burnin = 0, thin = 10)
chain3_beta1_thin<-thin_param_MCMC(chain3_beta1,burnin = 0, thin = 10)
chain3_psi_thin<-thin_param_MCMC(chain3_psi,burnin = 0, thin = 10)
chain3_tau_thin<-thin_param_MCMC(chain3_tau,burnin = 0, thin = 10)

# saving parallel chains
Beta0<- list(chain1_beta0_thin,chain2_beta0_thin,chain3_beta0_thin)
Beta1<- list(chain1_beta1_thin,chain2_beta1_thin,chain3_beta1_thin)
psi<-list(chain1_psi_thin,chain2_psi_thin,chain3_psi_thin)
tau<- list(chain1_tau_thin,chain2_tau_thin,chain3_tau_thin)

# checking convergence
library(coda)
library(mcmcplots)
Beta0.mcmc<-convert.mcmc.list(chain1_beta0)
Beta1.mcmc<-convert.mcmc.list(chain1_beta1)
Beta2.mcmc<-convert.mcmc.list(chain1_beta2)
tau.mcmc<-convert.mcmc.list(tau)
psi.mcmc<-convert.mcmc.list(psi)
g1<-gelman.diag(Beta0.mcmc)
g2<-gelman.diag(Beta1.mcmc)
g3<-gelman.diag(tau.mcmc)
g4<-gelman.diag(psi.mcmc)

# combining all chains to calculate summary and densities
Beta0_all<- c(chain1_beta0_thin,chain2_beta0_thin,chain3_beta0_thin)
Beta1_all<-c(chain1_beta1_thin,chain2_beta1_thin,chain3_beta1_thin)
tau_all<- c(chain1_tau_thin,chain2_tau_thin,chain3_tau_thin)
psi_all<-list()
for(i in 1:2148){
  psi_i<- c(chain1_psi_thin[i,],chain2_psi_thin[i,],chain3_psi_thin[i,])
  psi_all[[i]]<- psi_i
}
psi1_all<-c(chain1_psi_thin,chain2_psi_thin,chain3_psi_thin)
psi56_all<-c(chain1_psi_thin,chain2_psi_thin,chain3_psi_thin)

#saving results for later use
BYM_all_scotlip<-list(Beta0_all,Beta1_all,tau_all,psi_all)
BYM_mcmc_scotlip<-list(Beta0.mcmc,Beta1.mcmc,tau.mcmc,psi.mcmc)
save(BYM_all_scotlip,file="ind_all_ospg_male_surv.RData")
save(BYM_mcmc_scotlip,file="ind_mcmc_ospg_male_surv.RData")