source("Data preparation.R")
library("dplyr")
library("broom")
library(tidyverse)
head(dat)


# Axis limits and colours
getRange = function(x){
  vals = log(x) ###Confused where vals came from
  vals = vals[!is.na(vals)&is.finite(vals)]
  return(range(vals))
}

##Conor's code to work on generating plots editing in progress  
mitochan = "VDAC1"
mitorange = getRange(dat[,mitochan]) ####struggling with mitorange since I am unsure how to convert for mine

storeDefect = list()

pdf("Conf_Interval_95_Exp1.pdf",width=11.69,height=8.27)
# Do everything in a loop: always avoid copy-pasting code where you can: minimises mistakes and inconsistencies
for(pat in pats){
  datapat = subset(dat, patient_id==pat)
  allchans = gsub("LOG_","",colnames(dat)[grepl("LOG",colnames(dat))])
  chans = allchans
  chanvec = chans[chans%in%colnames(datapat)]
  chanvec = chanvec[chanvec!=mitochan]
  chandefect = rep(-999,length(chanvec))
  names(chandefect) = chanvec
  op = par(mfrow=c(3,4),mar=c(4.5,4.5,1,1))
  for(chan in chanvec){
    dt = subset(dat, patient_id==pat)
    chanrange = getRange(dt[,chan])
    # Linear regression for controls and predictive interval
    xctrl = log(datac[,mitochan])
    yctrl = log(datac[,chan])
    xpat = log(dt[,mitochan])
    ypat = log(dt[,chan])
    xsyn = seq(mitorange[1],mitorange[2],length.out=50)
    mod = lm(yctrl~xctrl)
  
    pred_syn = predict(mod,newdata = data.frame(xctrl=xsyn), se.fit=TRUE,  interval = "prediction",na.action=na.omit)$fit
    mid_syn = pred_syn[,1]
    up_syn = pred_syn[,3]
    low_syn = pred_syn[,2]
  
    pred = predict(mod,newdata = data.frame(xctrl=xpat), se.fit=TRUE,  interval = "prediction",na.action=na.omit)$fit
    up = pred[,3]
    low = pred[,2]	 
    dt[[paste(chan,"up_pred",sep="_")]]=up
    dt[[paste(chan,"low_pred",sep="_")]]=low
    dt[[paste(chan,"under_exp",sep="_")]] = log(dt[[chan]]) < dt[[paste(chan,"low_pred",sep="_")]]
    dt[[paste(chan,"over_exp",sep="_")]] = log(dt[[chan]]) > dt[[paste(chan,"up_pred",sep="_")]]
    
    ###Paste to dat
    datapat[[paste(chan,"up_pred",sep="_")]]=up
    datapat[[paste(chan,"low_pred",sep="_")]]=low
    datapat[[paste(chan,"under_exp",sep="_")]] = log(dt[[chan]]) < dt[[paste(chan,"low_pred",sep="_")]]
    datapat[[paste(chan,"over_exp",sep="_")]] = log(dt[[chan]]) > dt[[paste(chan,"up_pred",sep="_")]]
  
    #psd = (up - low)/(2*1.96)
    #upz = mid+3*psd
    #lowz = mid-3*psd 
    
    N = length(dt$ID)
    fracUP = sum(dt[[paste(chan,"over_exp",sep="_")]],na.rm=TRUE)/N
    fracDOWN = sum(dt[[paste(chan,"under_exp",sep="_")]],na.rm=TRUE)/N
    chandefect[chan] = fracDOWN
    
    mlab = paste(pat,"below:",signif(100*fracDOWN,3),"% above:",signif(100*fracUP,3),"%")
  
    plot(xpat,ypat,xlab=paste("log(",mitochan,")",sep=""),ylab=paste("log(",chan,")",sep=""),xlim=mitorange,ylim=chanrange,col=rgb(1,0,0,0.2),pch=16,main=mlab,cex.lab=2,cex.axis=1.75, abline(v=1.05, col="black", lty=4, lwd=1))
    points(xsyn,mid_syn,type="l",lty=1,lwd=2,col="black")
    points(xsyn,up_syn,type="l",lty=2,lwd=2,col="black")
    points(xsyn,low_syn,type="l",lty=2,lwd=2,col="black")
  }
  storeDefect[[pat]] = chandefect
  
  par(op)
mtext(paste("N =",N), outer=TRUE,  cex=1, line=-1.3)
write.csv(datapat, file=paste(pat,"_95PI_sig_data.csv",sep=""),row.names=FALSE)
}
dev.off()

defectdf = do.call("rbind",storeDefect)


print(storeDefect)

write.table(storeDefect, file="Summary95PI.csv",row.names=FALSE,sep=",")


####Bind all 95%PI classifications
P01dat = read.delim("P01_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
P02dat = read.delim("P02_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
P03dat = read.delim("P03_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
P04dat = read.delim("P04_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
P05dat = read.delim("P05_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
P06dat = read.delim("P06_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
P07dat = read.delim("P07_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
P09dat = read.delim("P09_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
P10dat = read.delim("P10_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
P11dat = read.delim("P11_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
P12dat = read.delim("P12_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
C01dat = read.delim("C01_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
C02dat = read.delim("C02_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
C03dat = read.delim("C03_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)
C04dat = read.delim("C04_95PI_sig_data.csv",sep=",",stringsAsFactors=FALSE)

Allpdat95PI <- rbind(P12dat,P06dat,C01dat,P11dat,C04dat,P10dat,P07dat,P09dat,P03dat,P04dat,P02dat,P01dat,C02dat)
write.table(Allpdat95PI, file="Allpdat95PI.csv",row.names=FALSE,sep=",")

####PLotting data by 95% predictive interval classification

Pred_Int<-read.csv("P01_95PI_sig_data.csv", header=T)
colnames(Pred_Int)

#Function SEM
sem = function(x){ 
  sd(x)/sqrt(length(x))
}
p = Pred_Int %>% select(5:62)
chans = colnames(p)

pdf("P01_isolatedCI_comp.pdf",width=6,height=8)
df=CIisol
op = par(mfrow=c(2,3))
for(chan in chans){
  Num = length(df$ID)
  my_mean = mean(df[[chan]])
  my_sem = sem(df[[chan]])
  stripchart(df[[chan]] ~ df$NDUFB8_under_exp,vertical=TRUE,method="jitter",jitter=0.1,pch=16,col=rgb(1,0,0,0.1),cex=1.25,
             main = paste("N =", Num), xlab = "Cluster", ylab = chan,cex.lab=1.35)
  boxplot(df[[chan]] ~ df$NDUFB8_under_exp, data = df,add=TRUE,outline=FALSE,notch=TRUE,col=NULL,boxwex=0.3,lwd=2)
}
par(op)
dev.off()

##Looking at isolated deficiencies
###Define deficiency groups NDUFB8 
print(chans)

#Filter for RRF
#CIisol= Pred_Int %>% filter(LOG_VDAC1>1.34)

#CIisol= Pred_Int %>% filter(MTCYB_under_exp==FALSE)
#CIisol= CIisol %>% filter(NDUFB8_under_exp==FALSE)
#CIisol= CIisol %>% filter(MTCO1_under_exp==FALSE)
#CIisol= CIisol %>% filter(SDHA_under_exp==FALSE)
#CIisol= CIisol %>% filter(NDUFB8_under_exp==TRUE)
#CIisol= CIisol %>% filter(ATP58_under_exp==FALSE)

#CIisol= CIisol %>% filter(MTND4_under_exp==FALSE)
#CIisol= CIisol %>% filter(UqCRC2_under_exp==FALSE)
#CIisol= CIisol %>% filter(COX4_under_exp==FALSE)
print(length(Pred_Int$ID))
print(length(CIisol$ID))


