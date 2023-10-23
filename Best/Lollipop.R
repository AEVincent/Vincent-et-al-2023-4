library("hrbrthemes")
library("dplyr")
library("ggplot2")

##Open all pats and bind together

P01 = read.delim("P01_catTri1.csv",sep=",",stringsAsFactors=FALSE)
P02 = read.delim("P02_catTri1.csv",sep=",",stringsAsFactors=FALSE)
P03 = read.delim("P03_catTri1.csv",sep=",",stringsAsFactors=FALSE)
P04 = read.delim("P04_catTri1.csv",sep=",",stringsAsFactors=FALSE)
P05 = read.delim("P05_catTri1.csv",sep=",",stringsAsFactors=FALSE)
P06 = read.delim("P06_catTri1.csv",sep=",",stringsAsFactors=FALSE)
P07 = read.delim("P07_catTri1.csv",sep=",",stringsAsFactors=FALSE)
P08 = read.delim("P08_catTri1.csv",sep=",",stringsAsFactors=FALSE)
P09 = read.delim("P09_catTri1.csv",sep=",",stringsAsFactors=FALSE)
P10 = read.delim("P10_catTri1.csv",sep=",",stringsAsFactors=FALSE)
P11 = read.delim("P11_catTri1.csv",sep=",",stringsAsFactors=FALSE)
P12 = read.delim("P12_catTri1.csv",sep=",",stringsAsFactors=FALSE)

All<-rbind(P01, P02, P03, P05, P06, P09, P10)

write.csv(All, file="AllPats_NormalvRRF.csv",row.names=FALSE)

###Read and plot
All = read.delim("AllPats_NonevCI+CIV.csv",sep=",",stringsAsFactors=FALSE)

cdPalette=c("darkorange","dodgerblue3","darkred","grey80")
level_order=c("None","CI+CIV")

#a=read.csv("asyn.csv")
head(All)

fig1=ggplot(All,aes(x=cat, y=muDiff, color=casecode))+
  geom_point(aes(size=abs(effSz)), alpha=0.5)+
  theme_bw() + 
  ggtitle("")+
  coord_flip()+  
  theme(
    axis.text=element_text(size=12),axis.title=element_text(size=16),
    plot.title = element_text(hjust = 0),
    panel.border = element_rect(colour = "black", fill=NA))+
  ylab("u1-u2_VDAC1")+
  xlab("")+  
  geom_hline(yintercept=0, linetype="dashed")
fig1

gd1 <- All%>%
  group_by(cat)%>% 
  mutate(median.value = mean(muDiff))

gd2 <- gd1%>%
  group_by(cat)%>% 
  mutate(sd = sd(muDiff))

fig1=ggplot(gd2,aes(x=reorder(cat,median.value),y=muDiff,size=abs(effSz)))+
  geom_point(alpha=0.8,aes(colour=factor(casecode),size=abs(effSz)))+
  ylim(-20, 20)+
  geom_segment(aes(x=reorder(cat,median.value), xend=reorder(cat,median.value), y=0, yend=median.value),size=3, alpha=0.03)+
  geom_errorbar(aes(ymin =median.value-sd, ymax = median.value+sd),size=0.5, width = 0.2,alpha=0.1)+
  theme_bw()+
  coord_flip()+  
  geom_hline(yintercept=0, linetype="dashed",colour="darkred")+
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=12))+
  ylab(expression(paste(mu, "1-", mu, "2")))+
  xlab("")+labs(size="Effect Size",col="")+
  ggtitle("CIvCIV")
fig1

jpeg("CIvCIVrank.jpg",width=450,height=500)
print(fig1)
dev.off()


