qpcr<-read.table("~/Projects/PHD_Host_Response_BV_Dynamics/Subprojects/PHD_miRNA/Experimental_Results/Migration/092316/092616_LCM_QPCR.txt",sep="\t",header=T)
library(ggplot2)
library(dplyr)
library(reshape)
ggplot(filter(qpcr,InputRNA==7 & InputcDNA==3))+geom_boxplot(aes(x=Sample,y=-Ct,col=PrimerPair))
ggplot(filter(qpcr,InputRNA==.7 | InputcDNA==.3))+geom_boxplot(aes(x=Sample,y=-Ct,col=PrimerPair))+facet_wrap(~InputRNA+InputcDNA)
scratch<-read.table("~/Projects/PHD_Host_Response_BV_Dynamics/Subprojects/PHD_miRNA/Experimental_Results/Migration/092316/AreaSummary.txt",sep="\t",header=T)
scratch<-melt(scratch)
ggplot(scratch)+geom_boxplot(aes(x=Sample,y=value))+theme_bw()+ggtitle("Percent VK2 Cell Area Filled at 8 hours post Culture Media Exposure")+xlab("Culture Media")+ylab("Percent Area Filled")

phColorTbl.hmp
