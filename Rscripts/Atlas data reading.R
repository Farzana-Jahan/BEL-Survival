#reading the atlas data (published and downloaded)
.libPaths("c:/software/Rpackages")
library(tidyverse)
library(dplyr)
#library(plyr)
atlas_data<-read_csv("Data/EHR Downloadable Data.csv")
#calculating Standard deviation of  SIR

#calculating  of all necessary values to calculate sd of  sir
atlas_data<-atlas_data%>%
  mutate(p025=(p025),p50=(p50),p975=(p975))
summary(atlas_data)

#transforming variables to calculate sd( SIR)
#atlas_data<-atlas_data%>%
# mutate(sd=(p975-p50)/1.96)
atlas_data<-atlas_data%>%
  mutate(sd=(p975-p025)/(2*1.96))
#atlas_data$sd<-atlas_data$sd_2

summary(atlas_data)
summary(atlas_data$sd)
summary(atlas_data$sd_2)
# take sd_2 as the standard deviation

#subsetting data for different cancers
#Oesophegael Cancer
cancer.ospg<-atlas_data%>%
  filter(Cancer_code==11)
cancer.ospg.male<-cancer.ospg%>%
  filter(Sex_code==1)
cancer.ospg.female<-cancer.ospg%>%
  filter(Sex_code==2)
cancer.ospg.persons<-cancer.ospg%>%
  filter(Sex_code==3)
cancer.ospg.male<-cancer.ospg.male%>%
  rename(c("p50"="y.ospg.male")) %>%
  rename(c("sd"= "sd.ospg.male"))
cancer.ospg.female<-cancer.ospg.female%>%
  rename(c("p50"="y.ospg.female")) %>%
  rename(c("sd"= "sd.ospg.female"))
cancer.ospg.persons<-cancer.ospg.persons%>%
  rename(c("p50"="y.ospg.persons")) %>%
  rename(c("sd"= "sd.ospg.persons"))
ospg.model.data<-inner_join(cancer.ospg.male,cancer.ospg.female,by="SA2_code")
ospg.model.data<-inner_join(ospg.model.data,cancer.ospg.persons,by="SA2_code")
ospg.model.data<-ospg.model.data%>%
  select(SA2_code,y.ospg.male,sd.ospg.male,y.ospg.female,sd.ospg.female,y.ospg.persons,sd.ospg.persons)
summary(ospg.model.data)
write_csv(ospg.model.data,"Data/ospg.data.csv")

#stomach cancer

cancer.stomach<-atlas_data%>%
  filter(Cancer_code==12)
cancer.stomach.male<-cancer.stomach%>%
  filter(Sex_code==1)
cancer.stomach.female<-cancer.stomach%>%
  filter(Sex_code==2)
cancer.stomach.persons<-cancer.stomach%>%
  filter(Sex_code==3)
cancer.stomach.male<-cancer.stomach.male%>%
  rename(c("p50"="y.stomach.male","sd"="sd.stomach.male"))
cancer.stomach.female<-cancer.stomach.female%>%
  rename(c("p50"="y.stomach.female", "sd"="sd.stomach.female"))
cancer.stomach.persons<-cancer.stomach.persons%>%
  rename(c("p50"="y.stomach.persons", "sd"="sd.stomach.persons"))
stomach.model.data<-inner_join(cancer.stomach.male,cancer.stomach.female,by="SA2_code")
stomach.model.data<-inner_join(stomach.model.data,cancer.stomach.persons,by="SA2_code")
stomach.model.data<-stomach.model.data%>%
  select(SA2_code,y.stomach.male,sd.stomach.male,y.stomach.female,sd.stomach.female,y.stomach.persons,sd.stomach.persons)
summary(stomach.model.data)
write_csv(stomach.model.data,"stomach.data.csv")

#colorectal cancer

cancer.colorectal<-atlas_data%>%
  filter(Cancer_code==14)
cancer.colorectal.male<-cancer.colorectal%>%
  filter(Sex_code==1)
cancer.colorectal.female<-cancer.colorectal%>%
  filter(Sex_code==2)
cancer.colorectal.persons<-cancer.colorectal%>%
  filter(Sex_code==3)
cancer.colorectal.male<-cancer.colorectal.male%>%
  rename(c("p50"="y.colorectal.male", "sd"="sd.colorectal.male"))
cancer.colorectal.female<-cancer.colorectal.female%>%
  rename(c("p50"="y.colorectal.female", "sd"="sd.colorectal.female"))
cancer.colorectal.persons<-cancer.colorectal.persons%>%
  rename(c("p50"="y.colorectal.persons","sd"="sd.colorectal.persons"))
colorectal.model.data<-inner_join(cancer.colorectal.male,cancer.colorectal.female,by="SA2_code")
colorectal.model.data<-inner_join(colorectal.model.data,cancer.colorectal.persons,by="SA2_code")
colorectal.model.data<-colorectal.model.data%>%
  select(SA2_code,y.colorectal.male,sd.colorectal.male,y.colorectal.female,sd.colorectal.female,y.colorectal.persons,sd.colorectal.persons)
summary(colorectal.model.data)
write_csv(colorectal.model.data,"colorectal.data.csv")

#liver cancer

cancer.liver<-atlas_data%>%
  filter(Cancer_code==18)
cancer.liver.male<-cancer.liver%>%
  filter(Sex_code==1)
cancer.liver.female<-cancer.liver%>%
  filter(Sex_code==2)
cancer.liver.persons<-cancer.liver%>%
  filter(Sex_code==3)
cancer.liver.male<-cancer.liver.male%>%
  rename(c("p50"="y.liver.male", "sd"="sd.liver.male"))
cancer.liver.female<-cancer.liver.female%>%
  rename(c("p50"="y.liver.female","sd"=" sd.liver.female"))
cancer.liver.persons<-cancer.liver.persons%>%
  rename(c("p50"="y.liver.persons", "sd"="sd.liver.persons"))
liver.model.data<-inner_join(cancer.liver.male,cancer.liver.female,by="SA2_code")
liver.model.data<-inner_join(liver.model.data,cancer.liver.persons,by="SA2_code")
liver.model.data<-liver.model.data%>%
  select(SA2_code,y.liver.male,sd.liver.male,y.liver.female,sd.liver.female,y.liver.persons,sd.liver.persons)
summary(liver.model.data)
write_csv(liver.model.data,"liver.data.csv")

#pncr cancer

cancer.pncr<-atlas_data%>%
  filter(Cancer_code==20)
cancer.pncr.male<-cancer.pncr%>%
  filter(Sex_code==1)
cancer.pncr.female<-cancer.pncr%>%
  filter(Sex_code==2)
cancer.pncr.persons<-cancer.pncr%>%
  filter(Sex_code==3)
cancer.pncr.male<-cancer.pncr.male%>%
  rename(c("p50"="y.pncr.male", "sd"="sd.pncr.male"))
cancer.pncr.female<-cancer.pncr.female%>%
  rename(c("p50"="y.pncr.female","sd"="sd.pncr.female"))
cancer.pncr.persons<-cancer.pncr.persons%>%
  rename(c("p50"="y.pncr.persons", sd="sd.pncr.persons"))
pncr.model.data<-inner_join(cancer.pncr.male,cancer.pncr.female,by="SA2_code")
pncr.model.data<-inner_join(pncr.model.data,cancer.pncr.persons,by="SA2_code")
pncr.model.data<-pncr.model.data%>%
  select(SA2_code,y.pncr.male,sd.pncr.male,y.pncr.female,sd.pncr.female,y.pncr.persons,sd.pncr.persons)
summary(pncr.model.data)
write_csv(pncr.model.data,"pncr.data.csv")

#lung cancer

cancer.lung<-atlas_data%>%
  filter(Cancer_code==23)
cancer.lung.male<-cancer.lung%>%
  filter(Sex_code==1)
cancer.lung.female<-cancer.lung%>%
  filter(Sex_code==2)
cancer.lung.persons<-cancer.lung%>%
  filter(Sex_code==3)
cancer.lung.male<-cancer.lung.male%>%
  rename(c("p50"="y.lung.male","sd"="sd.lung.male"))
cancer.lung.female<-cancer.lung.female%>%
  rename(c("p50"="y.lung.female", "sd"="sd.lung.female"))
cancer.lung.persons<-cancer.lung.persons%>%
  rename(c("p50"="y.lung.persons","sd"="sd.lung.persons"))
lung.model.data<-inner_join(cancer.lung.male,cancer.lung.female,by="SA2_code")
lung.model.data<-inner_join(lung.model.data,cancer.lung.persons,by="SA2_code")
lung.model.data<-lung.model.data%>%
  select(SA2_code,y.lung.male,sd.lung.male,y.lung.female,sd.lung.female,y.lung.persons,sd.lung.persons)
summary(lung.model.data)
write_csv(lung.model.data,"lung.data.csv")

#melanoma cancer

cancer.melanoma<-atlas_data%>%
  filter(Cancer_code==27)
cancer.melanoma.male<-cancer.melanoma%>%
  filter(Sex_code==1)
cancer.melanoma.female<-cancer.melanoma%>%
  filter(Sex_code==2)
cancer.melanoma.persons<-cancer.melanoma%>%
  filter(Sex_code==3)
cancer.melanoma.male<-cancer.melanoma.male%>%
  rename(c("p50"="y.melanoma.male","sd"="sd.melanoma.male"))
cancer.melanoma.female<-cancer.melanoma.female%>%
  rename(c("p50"="y.melanoma.female", "sd"="sd.melanoma.female"))
cancer.melanoma.persons<-cancer.melanoma.persons%>%
  rename(c("p50"="y.melanoma.persons","sd"=" sd.melanoma.persons"))
melanoma.model.data<-inner_join(cancer.melanoma.male,cancer.melanoma.female,by="SA2_code")
melanoma.model.data<-inner_join(melanoma.model.data,cancer.melanoma.persons,by="SA2_code")
melanoma.model.data<-melanoma.model.data%>%
  select(SA2_code,y.melanoma.male,sd.melanoma.male,y.melanoma.female,sd.melanoma.female,y.melanoma.persons,sd.melanoma.persons)
summary(melanoma.model.data)
write_csv(melanoma.model.data,"melanoma.data.csv")

#breast cancer

cancer.breast<-atlas_data%>%
  filter(Cancer_code==33)

cancer.breast.female<-cancer.breast%>%
  filter(Sex_code==2)


cancer.breast.female<-cancer.breast.female%>%
  rename(c("p50"="y.breast.female","sd"="sd.breast.female"))

breast.model.data<-cancer.breast.female
breast.model.data<-breast.model.data%>%
  select(SA2_code,y.breast.female,sd.breast.female)
summary(breast.model.data)
write_csv(breast.model.data,"breast.data.csv")

#cervical cancer

cancer.cervical<-atlas_data%>%
  filter(Cancer_code==35)

cancer.cervical.female<-cancer.cervical%>%
  filter(Sex_code==2)
head(cervical.model.data)
cervical.model.data<-cancer.cervical.female

cervical.model.data<-rename(cervical.model.data,c("p50"="y.cervical.female", "sd"="sd.cervical.female"))


cervical.model.data<-cervical.model.data%>%
  select(SA2_code,y.cervical.female,sd.cervical.female)
summary(cervical.model.data)
write_csv(cervical.model.data,"cervical.data.csv")
#uterine cancer

cancer.uterine<-atlas_data%>%
  filter(Cancer_code==36)

cancer.uterine.female<-cancer.uterine%>%
  filter(Sex_code==2)


cancer.uterine.female<-cancer.uterine.female%>%
  rename(c("p50"="y.uterine.female", "sd"="sd.uterine.female"))

uterine.model.data<-cancer.uterine.female
uterine.model.data<-uterine.model.data%>%
  select(SA2_code,y.uterine.female,sd.uterine.female)
summary(uterine.model.data)
write_csv(uterine.model.data,"uterine.data.csv")

#ovarian cancer

cancer.ovarian<-atlas_data%>%
  filter(Cancer_code==37)

cancer.ovarian.female<-cancer.ovarian%>%
  filter(Sex_code==2)


cancer.ovarian.female<-cancer.ovarian.female%>%
  rename(c("p50"="y.ovarian.female", "sd"="sd.ovarian.female"))

ovarian.model.data<-cancer.ovarian.female
ovarian.model.data<-ovarian.model.data%>%
  select(SA2_code,y.ovarian.female,sd.ovarian.female)
summary(ovarian.model.data)
write_csv(ovarian.model.data,"ovarian.data.csv")

#prostate cancer

cancer.prostate<-atlas_data%>%
  filter(Cancer_code==39)

cancer.prostate.male<-cancer.prostate%>%
  filter(Sex_code==1)


cancer.prostate.male<-cancer.prostate.male%>%
  rename(c("p50"="y.prostate.male","sd"="sd.prostate.male"))

prostate.model.data<-cancer.prostate.male
prostate.model.data<-prostate.model.data%>%
  select(SA2_code,y.prostate.male,sd.prostate.male)
summary(prostate.model.data)
write_csv(prostate.model.data,"prostate.data.csv")

#kidney cancer

cancer.kidney<-atlas_data%>%
  filter(Cancer_code==42)
cancer.kidney.male<-cancer.kidney%>%
  filter(Sex_code==1)
cancer.kidney.female<-cancer.kidney%>%
  filter(Sex_code==2)
cancer.kidney.persons<-cancer.kidney%>%
  filter(Sex_code==3)
cancer.kidney.male<-cancer.kidney.male%>%
  rename(c("p50"="y.kidney.male","sd"="sd.kidney.male"))
cancer.kidney.female<-cancer.kidney.female%>%
  rename(c("p50"="y.kidney.female", "sd"="sd.kidney.female"))
cancer.kidney.persons<-cancer.kidney.persons%>%
  rename(c("p50"="y.kidney.persons", "sd"="sd.kidney.persons"))
kidney.model.data<-inner_join(cancer.kidney.male,cancer.kidney.female,by="SA2_code")
kidney.model.data<-inner_join(kidney.model.data,cancer.kidney.persons,by="SA2_code")
kidney.model.data<-kidney.model.data%>%
  select(SA2_code,y.kidney.male,sd.kidney.male,y.kidney.female,sd.kidney.female,y.kidney.persons,sd.kidney.persons)
summary(kidney.model.data)
write_csv(kidney.model.data,"kidney.data.csv")
#brain cancer

cancer.brain<-atlas_data%>%
  filter(Cancer_code==45)
cancer.brain.male<-cancer.brain%>%
  filter(Sex_code==1)
cancer.brain.female<-cancer.brain%>%
  filter(Sex_code==2)
cancer.brain.persons<-cancer.brain%>%
  filter(Sex_code==3)
cancer.brain.male<-cancer.brain.male%>%
  rename(c("p50"="y.brain.male","sd"="sd.brain.male"))
cancer.brain.female<-cancer.brain.female%>%
  rename(c("p50"="y.brain.female","sd"="sd.brain.female"))
cancer.brain.persons<-cancer.brain.persons%>%
  rename(c("p50"="y.brain.persons", "sd"="sd.brain.persons"))
brain.model.data<-inner_join(cancer.brain.male,cancer.brain.female,by="SA2_code")
brain.model.data<-inner_join(brain.model.data,cancer.brain.persons,by="SA2_code")
brain.model.data<-brain.model.data%>%
  select(SA2_code,y.brain.male,sd.brain.male,y.brain.female,sd.brain.female,y.brain.persons,sd.brain.persons)
summary(brain.model.data)
write_csv(brain.model.data,"brain.data.csv")

#Thyroid cancer

cancer.thyroid<-atlas_data%>%
  filter(Cancer_code==48)
cancer.thyroid.male<-cancer.thyroid%>%
  filter(Sex_code==1)
cancer.thyroid.female<-cancer.thyroid%>%
  filter(Sex_code==2)
cancer.thyroid.persons<-cancer.thyroid%>%
  filter(Sex_code==3)
cancer.thyroid.male<-cancer.thyroid.male%>%
  rename(c("p50"="y.thyroid.male", "sd"="sd.thyroid.male"))
cancer.thyroid.female<-cancer.thyroid.female%>%
  rename(c("p50"="y.thyroid.female", "sd"="sd.thyroid.female"))
cancer.thyroid.persons<-cancer.thyroid.persons%>%
  rename(c("p50"="y.thyroid.persons","sd"="sd.thyroid.persons"))
thyroid.model.data<-inner_join(cancer.thyroid.male,cancer.thyroid.female,by="SA2_code")
thyroid.model.data<-inner_join(thyroid.model.data,cancer.thyroid.persons,by="SA2_code")
thyroid.model.data<-thyroid.model.data%>%
  select(SA2_code,y.thyroid.male,sd.thyroid.male,y.thyroid.female,sd.thyroid.female,y.thyroid.persons,sd.thyroid.persons)
summary(thyroid.model.data)
write_csv(thyroid.model.data,"thyroid.data.csv")

#nh_lymph cancer

cancer.nh_lymph<-atlas_data%>%
  filter(Cancer_code==53)
cancer.nh_lymph.male<-cancer.nh_lymph%>%
  filter(Sex_code==1)
cancer.nh_lymph.female<-cancer.nh_lymph%>%
  filter(Sex_code==2)
cancer.nh_lymph.persons<-cancer.nh_lymph%>%
  filter(Sex_code==3)
cancer.nh_lymph.male<-cancer.nh_lymph.male%>%
  rename(c("p50"="y.nh_lymph.male","sd"="sd.nh_lymph.male"))
cancer.nh_lymph.female<-cancer.nh_lymph.female%>%
  rename(c("p50"="y.nh_lymph.female", "sd"="sd.nh_lymph.female"))
cancer.nh_lymph.persons<-cancer.nh_lymph.persons%>%
  rename(c("p50"="y.nh_lymph.persons", "sd"="sd.nh_lymph.persons"))
nh_lymph.model.data<-inner_join(cancer.nh_lymph.male,cancer.nh_lymph.female,by="SA2_code")
nh_lymph.model.data<-inner_join(nh_lymph.model.data,cancer.nh_lymph.persons,by="SA2_code")
nh_lymph.model.data<-nh_lymph.model.data%>%
  select(SA2_code,y.nh_lymph.male,sd.nh_lymph.male,y.nh_lymph.female,sd.nh_lymph.female,y.nh_lymph.persons,sd.nh_lymph.persons)
summary(nh_lymph.model.data)
write_csv(nh_lymph.model.data,"nh_lymph.data.csv")

#leukaemia cancer

cancer.leukaemia<-atlas_data%>%
  filter(Cancer_code==54)
cancer.leukaemia.male<-cancer.leukaemia%>%
  filter(Sex_code==1)
cancer.leukaemia.female<-cancer.leukaemia%>%
  filter(Sex_code==2)
cancer.leukaemia.persons<-cancer.leukaemia%>%
  filter(Sex_code==3)
cancer.leukaemia.male<-cancer.leukaemia.male%>%
  rename(c("p50"="y.leukaemia.male","sd"="sd.leukaemia.male"))
cancer.leukaemia.female<-cancer.leukaemia.female%>%
  rename(c("p50"="y.leukaemia.female","sd"="sd.leukaemia.female"))
cancer.leukaemia.persons<-cancer.leukaemia.persons%>%
  rename(c("p50"="y.leukaemia.persons","sd"="sd.leukaemia.persons"))
leukaemia.model.data<-inner_join(cancer.leukaemia.male,cancer.leukaemia.female,by="SA2_code")
leukaemia.model.data<-inner_join(leukaemia.model.data,cancer.leukaemia.persons,by="SA2_code")
leukaemia.model.data<-leukaemia.model.data%>%
  select(SA2_code,y.leukaemia.male,sd.leukaemia.male,y.leukaemia.female,sd.leukaemia.female,y.leukaemia.persons,sd.leukaemia.persons)
summary(leukaemia.model.data)
write_csv(leukaemia.model.data,"leukaemia.data.csv")

#myeloma cancer

cancer.myeloma<-atlas_data%>%
  filter(Cancer_code==60)
cancer.myeloma.male<-cancer.myeloma%>%
  filter(Sex_code==1)
cancer.myeloma.female<-cancer.myeloma%>%
  filter(Sex_code==2)
cancer.myeloma.persons<-cancer.myeloma%>%
  filter(Sex_code==3)
cancer.myeloma.male<-cancer.myeloma.male%>%
  rename(c("p50"="y.myeloma.male","sd"="sd.myeloma.male"))
cancer.myeloma.female<-cancer.myeloma.female%>%
  rename(c("p50"="y.myeloma.female","sd"="sd.myeloma.female"))
cancer.myeloma.persons<-cancer.myeloma.persons%>%
  rename(c("p50"="y.myeloma.persons", "sd"="sd.myeloma.persons"))
myeloma.model.data<-inner_join(cancer.myeloma.male,cancer.myeloma.female,by="SA2_code")
myeloma.model.data<-inner_join(myeloma.model.data,cancer.myeloma.persons,by="SA2_code")
myeloma.model.data<-myeloma.model.data%>%
  select(SA2_code,y.myeloma.male,sd.myeloma.male,y.myeloma.female,sd.myeloma.female,y.myeloma.persons,sd.myeloma.persons)
summary(myeloma.model.data)
write_csv(myeloma.model.data,"myeloma.data.csv")

#all cancer

cancer.all<-atlas_data%>%
  filter(Cancer_code==64)
cancer.all.male<-cancer.all%>%
  filter(Sex_code==1)
cancer.all.female<-cancer.all%>%
  filter(Sex_code==2)
cancer.all.persons<-cancer.all%>%
  filter(Sex_code==3)
cancer.all.male<-cancer.all.male%>%
  rename(c("p50"="y.all.male","sd"="sd.all.male"))
cancer.all.female<-cancer.all.female%>%
  rename(c("p50"="y.all.female", "sd"="sd.all.female"))
cancer.all.persons<-cancer.all.persons%>%
  rename(c("p50"="y.all.persons", "sd"="sd.all.persons"))
all.model.data<-inner_join(cancer.all.male,cancer.all.female,by="SA2_code")
all.model.data<-inner_join(all.model.data,cancer.all.persons,by="SA2_code")
all.model.data<-all.model.data%>%
  select(SA2_code,y.all.male,sd.all.male,y.all.female,sd.all.female,y.all.persons,sd.all.persons)
summary(all.model.data)
write_csv(all.model.data,"all.data.csv")

#head_neck cancer

cancer.head_neck<-atlas_data%>%
  filter(Cancer_code==71)
cancer.head_neck.male<-cancer.head_neck%>%
  filter(Sex_code==1)
cancer.head_neck.female<-cancer.head_neck%>%
  filter(Sex_code==2)
cancer.head_neck.persons<-cancer.head_neck%>%
  filter(Sex_code==3)
cancer.head_neck.male<-cancer.head_neck.male%>%
  rename(c("p50"="y.head_neck.male","sd"="sd.head_neck.male"))
cancer.head_neck.female<-cancer.head_neck.female%>%
  rename(c("p50"="y.head_neck.female","sd"="sd.head_neck.female"))
cancer.head_neck.persons<-cancer.head_neck.persons%>%
  rename(c("p50"="y.head_neck.persons","sd"="sd.head_neck.persons"))
head_neck.model.data<-inner_join(cancer.head_neck.male,cancer.head_neck.female,by="SA2_code")
head_neck.model.data<-inner_join(head_neck.model.data,cancer.head_neck.persons,by="SA2_code")
head_neck.model.data<-head_neck.model.data%>%
  select(SA2_code,y.head_neck.male,sd.head_neck.male,y.head_neck.female,sd.head_neck.female,y.head_neck.persons,sd.head_neck.persons)
summary(head_neck.model.data)
write_csv(head_neck.model.data,"head_neck.data.csv")

a<- ggplot(data=all.model.data,aes(sample=y.all.male))+stat_qq()+stat_qq_line()+labs(main="Checking normality of all cancer survival (male)")+theme_bw()

b<- ggplot(data=all.model.data,aes(sample=y.all.female))+stat_qq()+stat_qq_line()+labs(
                      main="Checking normality of all cancer survival (male)")+theme_bw()

c<- ggplot(data=brain.model.data,aes(sample=y.brain.male))+stat_qq()+stat_qq_line()+labs(main="Checking normality of brain cancer survival (male)")+theme_bw()
d<- ggplot(data=brain.model.data,aes(sample=y.brain.female))+stat_qq()+stat_qq_line()+labs(
  main="Checking normality of brain cancer survival (male)")+theme_bw()
library(ggpubr)
ggarrange(a,b,c,d,nrow=2,ncol=2)





