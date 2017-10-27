qpcr<-read.table("~/Desktop/092616_LCM_QPCR.txt",sep="\t",header=T)
library(ggplot2)
library(dplyr)
library(reshape)
ggplot(filter(qpcr,InputRNA==7 & InputcDNA==3))+geom_boxplot(aes(x=Sample,y=-Ct,col=PrimerPair))
ggplot(filter(qpcr,InputRNA==.7 | InputcDNA==.3))+geom_boxplot(aes(x=Sample,y=-Ct,col=PrimerPair))+facet_wrap(~InputRNA+InputcDNA)


sm<-read.table("~/Desktop/092316_ScratchMovement",sep="\t",header=T)


ggplot(melt(sm))+geom_boxplot(aes(x=Sample,y=value))+theme_bw()+ggtitle("Percent Scratch Area Filled after 8 hours CM Exposure")+xlab("LCM")+ylab("% Area Filled")
