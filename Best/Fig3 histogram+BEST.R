library(plyr)
library(dplyr)
library(ggplot2)
install.packages("BEST")
library("BEST")
?BEST

#upload data and catagorize
w = read.csv("Merged_GMM_RRF.csv") 
colnames(w)
chan=colnames(w)

hasLOG = grepl("LOG_",chan)
Log= chan[(hasLOG)]

title <- colnames(w)[c(2,44)]

#Subset dataset
LOG=w%>%select(title,Log)
colnames(LOG)

#Renames
#ZS1=ZS[,c(3-27)]
#colnames(ZS1)
#level_order <- c("NDUFA13","SDHA","MTCOI","ATP5B","OSCP","Porin","TFAM","PGC1","LnP",
#                 "Prohibitin","Htra2","GPS2","DJ1","Ser65")
#names(ZS1)[2:15]=level_order
#colnames(ZS1)

#Subset for BEST
z=LOG%>%filter(SectionID=="P12")#mannually change
colnames(z)

z1= z%>% filter(Manual_RRF=="TRUE")# check the number of cells 3 or more
z2= z%>% filter(Manual_RRF=="FALSE")


#BEST analysis
Sys.time()

pdf("P12_CYB.pdf") 
y1=z1 %>% pull(LOG_CYB)
y2=z2 %>% pull(LOG_CYB)
BESTout <- BESTmcmc(y1, y2)# BEST modelling
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
LOG_CYB<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
                ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_VDAC.pdf") 
y2=z1 %>% pull(LOG_VDAC)
y1=z2 %>% pull(LOG_VDAC)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
LOG_VDAC<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_SDHA.pdf") 
y2=z1 %>% pull(LOG_SDHA)
y1=z2 %>% pull(LOG_SDHA)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
LOG_SDHA<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
                ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_MTCO1pdf") 
y2=z1 %>% pull(LOG_MTCO1)
y1=z2 %>% pull(LOG_MTCO1)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
LOG_MTCO1<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_NDUFB8.pdf") 
y2=z1 %>% pull(LOG_NDUFB8)
y1=z2 %>% pull(LOG_NDUFB8)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
LOG_NDUFB8<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
                ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_OSCP.pdf") 
y2=z1 %>% pull(LOG_OSCP)
y1=z2 %>% pull(LOG_OSCP)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
LOG_OSCP<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
                ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_ND4.pdf") 
y2=z1 %>% pull(LOG_ND4)
y1=z2 %>% pull(LOG_ND4)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
LOG_ND4<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
               ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_COX4.pdf") 
y2=z1 %>% pull(LOG_COX4)
y1=z2 %>% pull(LOG_COX4)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
LOG_COX4<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
                ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_ATP5B.pdf") 
y2=z1 %>% pull(LOG_ATP5B)
y1=z2 %>% pull(LOG_ATP5B)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
LOG_ATP5B<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_UqCRC2.pdf") 
y2=z1 %>% pull(LOG_UqCRC2)
y1=z2 %>% pull(LOG_UqCRC2)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
LOG_UqCRC2<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_Dyrstrophin.pdf") 
y2=z1 %>% pull(LOG_Dystrophin)
y1=z2 %>% pull(LOG_Dystrophin)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
LOG_Dystrophin<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()


Sys.time()

###Bind together

P12<-rbind(LOG_CYB,LOG_VDAC,LOG_SDHA,LOG_MTCO1,LOG_NDUFB8,LOG_OSCP,LOG_ND4,LOG_COX4,LOG_ATP5B,
           LOG_UqCRC2, LOG_Dystrophin)

names=rownames(P12)
P12s=as.data.frame(P12[,1])  
P12s=cbind(names,P12s)
P12s$casecode=c("P12")
P12s=P12s%>%filter(names=="muDiff"|names=="effSz")

write.table(P12s,file="P12s_catTri.csv",sep=",",row.names = F)
write.table(P12,file="P12_catTri.csv",sep=",")

