library(plyr)
library(dplyr)
library(ggplot2)
install.packages("BEST")
library("BEST")
?BEST

#upload data and catagorize
w = read.csv("GMM_Classif_SingPix.csv") 
colnames(w)
chan=colnames(w)

hasLOG = grepl("Nuc_",chan)
Log= chan[(hasLOG)]

title <- colnames(w)[c(2,87)]

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

z1= z%>% filter(Group=="None")# check the number of cells 3 or more
z2= z%>% filter(Group=="CI+CIV")


#BEST analysis
Sys.time()

pdf("P12_CYB.pdf") 
y1=z1 %>% pull(Nuc_CYB_mean)
y2=z2 %>% pull(Nuc_CYB_mean)
BESTout <- BESTmcmc(y1, y2)# BEST modelling
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_CYB_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
                ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_HexokinaseII.pdf") 
y2=z1 %>% pull(Nuc_HexokinaseII_mean)
y1=z2 %>% pull(Nuc_HexokinaseII_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_HexokinaseII_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_VDAC.pdf") 
y2=z1 %>% pull(Nuc_VDAC_mean)
y1=z2 %>% pull(Nuc_VDAC_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_VDAC_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
                ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_DJ1pdf") 
y2=z1 %>% pull(Nuc_DJ1_mean)
y1=z2 %>% pull(Nuc_DJ1_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_DJ1_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_PDH.E2.pdf") 
y2=z1 %>% pull(Nuc_PDH.E2_mean)
y1=z2 %>% pull(Nuc_PDH.E2_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_PDH.E2_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
                ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_HSP60.pdf") 
y2=z1 %>% pull(Nuc_HSP60_mean)
y1=z2 %>% pull(Nuc_HSP60_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_HSP60_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
                ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_MTHFD2.pdf") 
y2=z1 %>% pull(Nuc_MTHFD2_mean)
y1=z2 %>% pull(Nuc_MTHFD2_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_MTHFD2_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
               ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_GPS2.pdf") 
y2=z1 %>% pull(Nuc_GPS2_mean)
y1=z2 %>% pull(Nuc_GPS2_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_GPS2_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
                ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_TFAM.pdf") 
y2=z1 %>% pull(Nuc_TFAM_mean)
y1=z2 %>% pull(Nuc_TFAM_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_TFAM_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_MTCO1.pdf") 
y2=z1 %>% pull(Nuc_MTCO1_mean)
y1=z2 %>% pull(Nuc_MTCO1_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_MTCO1_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_SIRT3.pdf") 
y2=z1 %>% pull(Nuc_SIRT3_mean)
y1=z2 %>% pull(Nuc_SIRT3_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_SIRT3_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_NDUFB8.pdf") 
y2=z1 %>% pull(Nuc_NDUFB8_mean)
y1=z2 %>% pull(Nuc_NDUFB8_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_NDUFB8_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_ATF4.pdf") 
y2=z1 %>% pull(Nuc_ATF4_mean)
y1=z2 %>% pull(Nuc_ATF4_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_ATF4_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_ATF6.pdf") 
y2=z1 %>% pull(Nuc_ATF6_mean)
y1=z2 %>% pull(Nuc_ATF6_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_ATF6_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_Htra2.pdf") 
y2=z1 %>% pull(Nuc_Htra2_mean)
y1=z2 %>% pull(Nuc_Htra2_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_Htra2_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_ClpP.pdf") 
y2=z1 %>% pull(Nuc_ClpP_mean)
y1=z2 %>% pull(Nuc_ClpP_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_ClpP_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_CHOP.pdf") 
y2=z1 %>% pull(Nuc_CHOP_mean)
y1=z2 %>% pull(Nuc_CHOP_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_CHOP_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()


pdf("P12_DRP1.pdf") 
y2=z1 %>% pull(Nuc_DRP1_mean)
y1=z2 %>% pull(Nuc_DRP1_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_DRP1_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_Beclin1.pdf") 
y2=z1 %>% pull(Nuc_Beclin1_mean)
y1=z2 %>% pull(Nuc_Beclin1_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_Beclin1_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_Prohibitin.pdf") 
y2=z1 %>% pull(Nuc_Prohibitin_mean)
y1=z2 %>% pull(Nuc_Prohibitin_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_Prohibitin_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

pdf("P12_LONP1.pdf") 
y2=z1 %>% pull(Nuc_LONP1_mean)
y1=z2 %>% pull(Nuc_LONP1_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_LONP1_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()


pdf("P12_PGC1a.pdf") 
y2=z1 %>% pull(Nuc_PGC1a_mean)
y1=z2 %>% pull(Nuc_PGC1a_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_PGC1a_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()


pdf("P12_Dystrophin.pdf") 
y2=z1 %>% pull(Nuc_Dystrophin_mean)
y1=z2 %>% pull(Nuc_Dystrophin_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_Dystrophin_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()


pdf("P12_DNA1.pdf") 
y2=z1 %>% pull(Nuc_DNA1_mean)
y1=z2 %>% pull(Nuc_DNA1_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_DNA1_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()


pdf("P12_DNA2.pdf") 
y2=z1 %>% pull(Nuc_DNA2_mean)
y1=z2 %>% pull(Nuc_DNA2_mean)
BESTout <- BESTmcmc(y1, y2)
plotAll(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), 
        ROPEeff=c(-0.2,0.2), compValm=0.5, showCurve=TRUE) 
Nuc_DNA2_mean<-summary(BESTout, credMass=0.95, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
              ROPEeff=c(-0.2,0.2))
dev.off()

Sys.time()

###Bind together

P12<-rbind(Nuc_CYB_mean,Nuc_HexokinaseII_mean,Nuc_VDAC_mean,Nuc_DJ1_mean,Nuc_PDH.E2_mean,Nuc_HSP60_mean,Nuc_MTHFD2_mean,Nuc_GPS2_mean,Nuc_TFAM_mean,Nuc_MTCO1_mean,
           Nuc_SIRT3_mean,Nuc_NDUFB8_mean,Nuc_ATF4_mean,Nuc_ATF6_mean,Nuc_Htra2_mean,Nuc_ClpP_mean,Nuc_CHOP_mean,Nuc_DRP1_mean,Nuc_Beclin1_mean,Nuc_Prohibitin_mean,Nuc_LONP1_mean,Nuc_PGC1a_mean,Nuc_Dystrophin_mean,Nuc_DNA1_mean,Nuc_DNA2_mean)

names=rownames(P12)
P12s=as.data.frame(P12[,1])  
P12s=cbind(names,P12s)
P12s$casecode=c("P12")
P12s=P12s%>%filter(names=="muDiff"|names=="effSz")

write.table(P12s,file="P12s_catTri.csv",sep=",",row.names = F)
write.table(P12,file="P12_catTri.csv",sep=",")

