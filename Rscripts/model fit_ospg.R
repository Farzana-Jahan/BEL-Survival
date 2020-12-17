# exploring data
#ospg survival male
ospg.model.data<- read_csv("Data/ospg.data.csv")
head(ospg.model.data)
library(tidyverse)
ggplot(data=ospg.model.data, aes(x= y.ospg.male))+geom_density()+theme_bw()

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
ospg.model.data<- inner_join(ospg.model.data,cov,by="SA2_code")
write_csv(ospg.model.data,"Data/ospg.data.csv")
# Parametric Model fit

library(CARBayes)
param_BYM<-S.CARleroux(y.ospg.male~region-1, data=ospg.model.data,W= W, family="gaussian",rho=1, n.sample = 100000,
                       burnin = 10000,thin=1)
save(param_BYM,"Results/ospg_BYM")
param_Ind<-S.CARleroux(y.ospg.male~region-1, data=ospg.model.data,W= W, family="gaussian",rho=0, n.sample = 100000,
                       burnin = 10000,thin=1)
save(param_Ind,"Results/ospg_Ind")
param_Leroux<-S.CARleroux(y.ospg.male~region-1, data=ospg.model.data,W= W, family="gaussian", n.sample = 100000,
                       burnin = 10000,thin=1)
save(param_Leroux,"Results/ospg_Leroux")
param_Leroux$modelfit
library(ngspatial)
param_moranBasis<-sparse.sglmm(y.ospg.male~region-1,data=ospg.model.data,A= W,
                               family="gaussian",minit = 10000, maxit=100000)
save(param_moranBasis,"Results/ospg_moranBasis")


