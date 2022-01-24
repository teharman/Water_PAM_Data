setwd("D:/PhD_PAM/Data_CSVs")

library(lme4)
library(sjPlot)
library(sjmisc)
library(glmmTMB)
library(ggplot2)
library(tidyverse)
library(reshape)
library(ggpubr)
library(Rmisc)

All_Data1<-read.csv("Main_Data.csv",header=T)

source("http://goo.gl/UUyEzD") #outlier code source
outlierKD(All_Data1, Yield) #outlier code

hist(All_Data1$Yield)

All_Data1$transplant.treatment=(paste(All_Data1$Translocation,All_Data1$Treatment,sep="."))
All_Data_E1<-subset(All_Data1,Day == "Experiment",select=c(Species,Fragment,transplant.treatment,Time,F,Fm.,Yield))
All_Data_E1<-subset(All_Data_E1,Species=="Porites astreoides",select=c(Fragment,transplant.treatment,Time,F,Fm.,Yield))
All_Data_EM1<-subset(All_Data_E1,Time == "Dusk",select=c(Fragment,transplant.treatment,F,Fm.,Yield))
All_Data_EM1<-All_Data_EM1[-c(6,28,54,58,70,119,120,121,171,194,222), ]
All_Data_EQ1<-subset(All_Data_E1,Time == "Noon",select=c(Fragment,transplant.treatment,F,Fm.,Yield))
All_Data_EQ1<-All_Data_EQ1[-c(13,15,34,70,121,143,147,152,153,188,251,253), ]

All_Data_E2<-subset(All_Data1,Day == "Experiment",select=c(Species,Fragment,Translocation,Treatment,Time,F,Fm.,Yield))
All_Data_E2<-subset(All_Data_E2,Species=="Porites astreoides",select=c(Fragment,Translocation,Treatment,Time,F,Fm.,Yield))
All_Data_EM2<-subset(All_Data_E2,Time == "Dusk",select=c(Fragment,Translocation,Treatment,F,Fm.,Yield))
All_Data_EM2<-All_Data_EM2[-c(6,28,54,58,70,119,120,121,171,194,222), ]
All_Data_EQ2<-subset(All_Data_E2,Time == "Noon",select=c(Fragment,Translocation,Treatment,F,Fm.,Yield))
All_Data_EQ2<-All_Data_EQ2[-c(13,15,34,70,121,143,147,152,153,188,251,253), ]

data.lm2<-glm(Yield~transplant.treatment+(1|Fragment),data=All_Data_EM1,family=inverse.gaussian(link="1/mu^2"))
summary(data.lm2) #no influence of treatment and translocation site on maximum quantum yield
plot(data.lm2)

data.lm3<-glm(Yield~transplant.treatment+(1|Fragment),data=All_Data_EQ1,family=inverse.gaussian(link="1/mu^2"))
summary(data.lm3) #influence of treatment and translocation site on effective quantum yield
plot(data.lm3)


#Now to look into relevant comparisons via lsmeans analysis

library(lsmeans)

PWC1<-lsmeans(data.lm2, specs = c("transplant.treatment"))
summary(PWC1)

rbind(pairs(PWC1), adjust = "none")

gd1<-summarySE(All_Data_EM2,measurevar="Yield", groupvars=c("Treatment","Translocation"))

PAMFigM <- ggplot(gd1, aes(x=Treatment, y=Yield, colour=Translocation,fill=Translocation,group=Translocation,shape=Translocation)) +
  geom_point(size=3,alpha=1,size = 3,stroke=7,position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5),size=1,linetype="longdash",alpha=1)+
  geom_errorbar(gd1,mapping=aes(ymin=Yield-ci, ymax=Yield+ci,color=Translocation,group=Translocation),position=position_dodge(width=0.5),size=1, alpha=1)+
  theme_bw(base_size=22)+
  theme(legend.spacing.y = unit(1.0, 'cm'),
        legend.key.size = unit(1.5, "cm"))+
  scale_x_discrete(breaks=c("28C","32C","32C+LPS"),
                   labels=c("28°C","32°C","32°C+LPS"))+
  labs(x="Treatment",y="Maximum Quantum Yield (Fv/Fm)")+
  theme(axis.title.y = element_text(size=20,vjust=4),
        axis.title.x = element_text(size=20,vjust=-2))+
  theme(axis.text.y = element_text(size=12))+
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
PAMFigM + geom_point(All_Data_EM2,mapping=aes(x=Treatment,y=Yield,colour=Translocation),alpha=0.2,shape=20,size=4,position = position_dodge(width = 0.5))


PWC2<-emmeans(data.lm3, specs = c("transplant.treatment"))
summary(PWC2)

rbind(pairs(PWC2), adjust = "none")

gd2<-summarySE(All_Data_EQ2,measurevar="Yield", groupvars=c("Treatment","Translocation"))

PAMFigQ <- ggplot(gd2, aes(x=Treatment, y=Yield, colour=Translocation,fill=Translocation,group=Translocation,shape=Translocation)) +
  geom_point(size=3,alpha=1,size = 5,stroke=7,position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5),size=1,linetype="longdash")+
  geom_errorbar(gd2,mapping=aes(ymin=Yield-ci, ymax=Yield+ci,color=Translocation,group=Translocation),position=position_dodge(width=0.5),size=1)+
  theme_bw(base_size=22)+
  theme(legend.spacing.y = unit(1.0, 'cm'),
        legend.key.size = unit(1.5, "cm"))+
  scale_x_discrete(breaks=c("28C","32C","32C+LPS"),
                   labels=c("28°C","32°C","32°C+LPS"))+
  labs(x="Treatment",y="Effective Quantum Yield (Fv'/Fm')")+
  theme(axis.title.y = element_text(size=20,vjust=4),
        axis.title.x = element_text(size=20,vjust=-2))+
  theme(axis.text.y = element_text(size=12))+
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
PAMFigQ + geom_point(All_Data_EQ2,mapping=aes(x=Treatment,y=Yield,colour=Translocation),alpha=0.2,shape=20,size=4,position = position_dodge(width = 0.5))
