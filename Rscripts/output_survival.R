# loading result of ind SBEL for ospg survival

load("U:/Research/Projects/sef/bayes_meta_el/SBEL_BYM_ospg_male_surv_1000.RData")
SBEL_ind_ospg_male_surv<-SBEL_BYM_ospg_male_surv
par(mfrow=c(2,2))
plot(1:1000,SBEL_ind_ospg_male_surv[[1]]$Beta[1,],type="l",col="green", main= "Traceplot for coeffcient(major cities)")
lines(1:1000,SBEL_ind_ospg_male_surv[[2]]$Beta[1,],type="l",col="blue")
lines(1:1000,SBEL_ind_ospg_male_surv[[3]]$Beta[1,],type="l",col="red")
plot(1:1000,SBEL_ind_ospg_male_surv[[1]]$Beta[2,],type="l",col="green", main= "Traceplot for coeffcient(regional)")
lines(1:1000,SBEL_ind_ospg_male_surv[[2]]$Beta[2,],type="l",col="blue")
lines(1:1000,SBEL_ind_ospg_male_surv[[3]]$Beta[2,],type="l",col="red")
plot(1:1000,SBEL_ind_ospg_male_surv[[1]]$Beta[3,],type="l",col="green", main= "Traceplot for coeffcient(remote)")
lines(1:1000,SBEL_ind_ospg_male_surv[[2]]$Beta[3,],type="l",col="blue")
lines(1:1000,SBEL_ind_ospg_male_surv[[3]]$Beta[3,],type="l",col="red")
plot(1:1000,SBEL_ind_ospg_male_surv[[1]]$tau,type="l",col="green", main= "Traceplot for precision")
lines(1:1000,SBEL_ind_ospg_male_surv[[2]]$tau,type="l",col="blue")
lines(1:1000,SBEL_ind_ospg_male_surv[[3]]$tau,type="l",col="red")

par(mfrow=c(2,2))
plot(density(SBEL_ind_ospg_male_surv[[1]]$Beta[1,]),col="green", main= "Posterior density for coeffcient(major cities)")
lines(density(SBEL_ind_ospg_male_surv[[2]]$Beta[1,]),col="blue")
lines(density(SBEL_ind_ospg_male_surv[[3]]$Beta[1,]),col="red")
plot(density(SBEL_ind_ospg_male_surv[[1]]$Beta[2,]),col="green", main= "Posterior density for coeffcient(regional)")
lines(density(SBEL_ind_ospg_male_surv[[2]]$Beta[2,]),col="blue")
lines(density(SBEL_ind_ospg_male_surv[[3]]$Beta[2,]),col="red")
plot(density(SBEL_ind_ospg_male_surv[[1]]$Beta[3,]),col="green", main= "Posterior density for coeffcient(remote)")
lines(density(SBEL_ind_ospg_male_surv[[2]]$Beta[3,]),col="blue")
lines(density(SBEL_ind_ospg_male_surv[[3]]$Beta[3,]),col="red")
plot(density(SBEL_ind_ospg_male_surv[[1]]$tau),col="green", main= "Posterior density for precision")
lines(density(SBEL_ind_ospg_male_surv[[2]]$tau),col="blue")
lines(density(SBEL_ind_ospg_male_surv[[3]]$tau),col="red")


# applying thining by 1 and burnin=1000 and saving as mcmc

thin_MCMC<- function(output,burnin=1000,thin=10)
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
    niter=dim(output)[2]
    out_postburn<-output[,(burnin+1):niter]
    n<-dim(out_postburn)[2]
    out_thin<- out_postburn[,seq(1,n,by=thin)]
    return(out_thin)
  }
  out_thin
}

chain1_beta0<-SBEL_ind_ospg_male_surv[[1]]$Beta[1,]
chain1_beta1<-SBEL_ind_ospg_male_surv[[1]]$Beta[2,]
chain1_beta2<-SBEL_ind_ospg_male_surv[[1]]$Beta[3,]
chain1_psi<-SBEL_ind_ospg_male_surv[[1]]$psi
chain1_tau<-SBEL_ind_ospg_male_surv[[1]]$tau
chain2_beta0<-SBEL_ind_ospg_male_surv[[2]]$Beta[1,]
chain2_beta1<-SBEL_ind_ospg_male_surv[[2]]$Beta[2,]
chain2_beta2<-SBEL_ind_ospg_male_surv[[2]]$Beta[3,]
chain2_psi<-SBEL_ind_ospg_male_surv[[2]]$psi
chain2_tau<-SBEL_ind_ospg_male_surv[[2]]$tau
chain3_beta0<-SBEL_ind_ospg_male_surv[[3]]$Beta[1,]
chain3_beta1<-SBEL_ind_ospg_male_surv[[3]]$Beta[2,]
chain3_beta2<-SBEL_ind_ospg_male_surv[[3]]$Beta[3,]
chain3_psi<-SBEL_ind_ospg_male_surv[[3]]$psi
chain3_tau<-SBEL_ind_ospg_male_surv[[3]]$tau

# Apply these to all BEL models after thinning and save with proper names
# thinned chains
chain1_beta0_thin<-thin_MCMC(chain1_beta0,burnin = 100, thin = 1)
chain1_beta1_thin<-thin_MCMC(chain1_beta1,burnin = 100, thin = 1)
chain1_beta2_thin<-thin_MCMC(chain1_beta2,burnin = 100, thin = 1)
chain1_psi_thin<-thin_MCMC(chain1_psi,burnin = 100, thin = 1)
chain1_tau_thin<-thin_MCMC(chain1_tau,burnin = 100, thin = 1)
chain2_beta0_thin<-thin_MCMC(chain2_beta0,burnin = 100, thin = 1)
chain2_beta1_thin<-thin_MCMC(chain2_beta1,burnin = 100, thin = 1)
chain2_beta2_thin<-thin_MCMC(chain2_beta2,burnin = 100, thin = 1)
chain2_psi_thin<-thin_MCMC(chain2_psi,burnin = 100, thin = 1)
chain2_tau_thin<-thin_MCMC(chain2_tau,burnin = 100, thin = 1)
chain3_beta0_thin<-thin_MCMC(chain3_beta0,burnin = 100, thin = 1)
chain3_beta1_thin<-thin_MCMC(chain3_beta1,burnin = 100, thin = 1)
chain3_beta2_thin<-thin_MCMC(chain3_beta2,burnin = 100, thin = 1)
chain3_psi_thin<-thin_MCMC(chain3_psi,burnin = 100, thin = 1)
chain3_tau_thin<-thin_MCMC(chain3_tau,burnin = 100, thin = 1)

# saving parallel chains
Beta0<- list(chain1_beta0_thin,chain2_beta0_thin,chain3_beta0_thin)
Beta1<- list(chain1_beta1_thin,chain2_beta1_thin,chain3_beta1_thin)
Beta2<- list(chain1_beta2_thin,chain2_beta2_thin,chain3_beta2_thin)
psi<-list(chain1_psi_thin,chain2_psi_thin,chain3_psi_thin)
tau<- list(chain1_tau_thin,chain2_tau_thin,chain3_tau_thin)

# checking convergence
library(coda)
library(mcmcplots)
Beta0.mcmc<-convert.mcmc.list(Beta0)
Beta1.mcmc<-convert.mcmc.list(Beta1)
Beta2.mcmc<-convert.mcmc.list(Beta2)
tau.mcmc<-convert.mcmc.list(tau)
psi.mcmc<-convert.mcmc.list(psi)
g1<-gelman.diag(Beta0.mcmc)
g2<-gelman.diag(Beta1.mcmc)
g3<-gelman.diag(Beta2.mcmc)
g4<-gelman.diag(tau.mcmc)
#g5<-gelman.diag(psi.mcmc)

# combining all chains to calculate summary and densities
Beta0_all<- c(chain1_beta0_thin,chain2_beta0_thin,chain3_beta0_thin)
Beta1_all<-c(chain1_beta1_thin,chain2_beta1_thin,chain3_beta1_thin)
Beta2_all<-c(chain1_beta2_thin,chain2_beta2_thin,chain3_beta2_thin)
tau_all<- c(chain1_tau_thin,chain2_tau_thin,chain3_tau_thin)
psi_all<-list()
for(i in 1:2148){
  psi_i<- c(chain1_psi_thin[i,],chain2_psi_thin[i,],chain3_psi_thin[i,])
  psi_all[[i]]<- psi_i
}
#psi1_all<-c(chain1_psi_thin,chain2_psi_thin,chain3_psi_thin)
#psi56_all<-c(chain1_psi_thin,chain2_psi_thin,chain3_psi_thin)

#saving results for later use
SBEL_ind_all_ospg_male_surv_male<-list(Beta0_all,Beta1_all,Beta2_all,tau_all,psi_all)
SBEL_ind_mcmc_ospg_male_surv<-list(Beta0.mcmc,Beta1.mcmc,Beta2.mcmc,tau.mcmc,psi.mcmc)
save(SBEL_ind_all_ospg_male_surv_male,file="Results/SBEL_ind_all_ospg_male_surv_male.RData")
save(SBEL_ind_mcmc_ospg_male_surv,file="Results/SBEL_ind_mcmc_ospg_male_surv.RData")

out<-data.frame(Paramter=c("Intercept", "Beta1(cities)", "Beta2(regional)","Precision"),
                Median=c(median(Beta0_all),median(Beta1_all),median(Beta2_all),median(tau_all)),
                `2.5% quantile`= c(quantile(Beta0_all,0.025),quantile(Beta1_all,0.025),quantile(Beta2_all,0.025),quantile(tau_all,0.025)),
                `97.5% quantile`= c(quantile(Beta0_all,0.975),quantile(Beta1_all,0.975),quantile(Beta2_all,0.975),quantile(tau_all,0.975)),
                `Gelman_Rubin diagnostic`= round(c(g1[[1]][1],g2[[1]][1],g3[[1]][1],g4[[1]][1]),2))
knitr::kable(out,row.names = F,format="latex")


# calculation of WAIC
source("Rscripts/waic.R")
Beta_ind_SBEL<- rbind(SBEL_ind_all_ospg_male_surv_male[[1]],SBEL_ind_all_ospg_male_surv_male[[2]],SBEL_ind_all_ospg_male_surv_male[[3]])
psi_ind_SBEL<-matrix(c(SBEL_ind_all_ospg_male_surv_male[[5]][[1]]),nrow=2700, ncol=1)
for(i in 2:2148){
  psi_ind_SBEL<-cbind(psi_ind_SBEL,SBEL_ind_all_ospg_male_surv_male[[4]][[i]])
}
tau_ind_BEL<-SBEL_ind_all_ospg_male_surv_male[[4]]
theta<-t(Beta_ind_SBEL)%*%t(x)+psi_ind_SBEL
#pred_BEL_ind<-colMeans(theta)
WAIC_ind_BEL<-get.WAIC.BEL(theta=theta,y=y,x=x)
