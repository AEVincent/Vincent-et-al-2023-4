source("Data preparation.R")
source("gmm_cluster.R")
library(MASS)

multiD_cluster = function(sec, Nclust, cluster_chans = c("MTCO1","VDAC1","NDUFB8")){
  
  sec_wide = reshape(sec[,c("Corrected","Chan","ID","SectionID")], #picks data from all row and selected channels
                     idvar = c("ID","SectionID"), #define identifiers of row, ID for fibres and SectionID for sections/regions/observer
                     timevar = "Chan", direction="wide")
  colnames(sec_wide) = gsub("Corrected.","",colnames(sec_wide))
  sec_w = log(sec_wide[,cluster_chans[cluster_chans%in%colnames(sec_wide)]])
  keepsec = !is.na(rowSums(sec_w))
  sec_w = sec_w[keepsec,]
  sec_wide = sec_wide[keepsec,]
  
  mb = Mclust(sec_w, Nclust)
  
  bc = mb$classification
  names(bc) = sec_wide$ID
  return(bc)
}

ctrl_ids = c("C01","C02","C03","C04")

pdat$Corrected = pdat$Value
pdat$Chan = pdat$Channel
pdat$SectionID = pdat$patient_id

##MITO DATASET
oxp = c("VDAC1","NDUFB8", "MTCYB", "MTCO1", "ATP5B", "COX4", "MTND4", "OSCP", "SDHA", "UqCRC2")
#oxp = c("VDAC1","MTND2", "MTCYB", "MTCO1", "ATP8", "COX4", "MTND4", "SDHA", "UqCRC2", "ATP5B")
CI = c("VDAC1", "NDUFB8", "MTND4")
NDUFB8 = c("VDAC1", "NDUFB8")
#CI = c("VDAC1", "MTND2", "MTND4")
#MTND2 = c("VDAC1", "MTND2")
MTND4 = c("VDAC1", "MTND4")
#CIII = c("VDAC1", "MTCYB")
CIII = c("VDAC1", "MTCYB", "UqCRC2")
MTCYB = c ("VDAC1", "MTCYB")
UQCRC2 = c("VDAC1", "UqCRC2")
CIV = c("VDAC1", "MTCO1", "COX4")
MTCO1 = c("VDAC1", "MTCO1")
COX4 = c("VDAC1", "COX4")
CV = c("VDAC1", "ATP5B", "OSCP")
#CV = c("VDAC1", "ATP8", "ATP5B")
OSCP = c("VDAC1", "OSCP")
ATP5B = c("VDAC1", "ATP5B")
#ATP8 = c("VDAC1", "ATP8")

##SIGNALLING DATASET
#oxp = c("VDAC1","NDUFB8", "MTCYB", "MTCO1")
#CI = c("VDAC1", "NDUFB8")
#CIII = c("VDAC1", "MTCYB")
#morph = c("Area","AspectRatio","Perimeter","Circularity","xCoord","yCoord")
#CIV = c("VDAC1", "MTCO1")
#nolog = chans[!grepl("LOG",chans)]
#justprots = nolog[!nolog%in%c("Mean","Edges",morph)]
#Sig1 = c("Htra2", "HexokinaseII")
#Sig2 = c("Hsp60", "LONP1", "PDHE2", "SIRT3", "Prohibitin", "TFAM", "ATF4", "ATF6", "Beclin")
#UPR = c("HSP60", "Htra2", "CHOP", "ClpP")
#Biogen = c("PGC1a", "SIRT3", "TFAM")
#Mitopro = c("Prohibitin", "TFAM", "LONP1")

ids = sort(unique(pdat$patient_id))
patids = ids[!ids%in%ctrl_ids]

#for (btypes in c("oxphos","cI", "cIII", "cIV", "signal","oxsig","oxsigmorph", "sig1", "sig2", "upr", "biogen", "mitopro")){
for (btypes in c("oxphos","cI","ndufb8", "mtnd4", "cIII", "mtcyb", "uqcrc2", "cIV", "mtco1", "cox4", "cV", "oscp", "atp5b")){
#for (btypes in c("oxphos","cI","MTND2", "MTND4", "cIII", "cIV", "mtCO1", "cox4", "cIV", "atp8", "atp5b")){
#for (btypes in c("oxphos")){
    
patlist = list()

for(patid in patids){
  print(btypes)
  print(patid)
  #patid = patids[1]
  pat = pdat[pdat$patient_id==patid,]
  ctrls = pdat[pdat$patient_id%in%ctrl_ids,]
  
  # OXPHOS only
  if (btypes=="oxphos") biochem = gmm_cluster(ctrls,pat,cluster_chans=oxp)
  # CI only
  if (btypes=="cI") biochem = gmm_cluster(ctrls,pat,cluster_chans=CI)
  # NDUFB8 only
  if (btypes=="ndufb8") biochem = gmm_cluster(ctrls,pat,cluster_chans=NDUFB8)
  # MTND2 only
  if (btypes=="mtnd2") biochem = gmm_cluster(ctrls,pat,cluster_chans=MTND2)
  # MTND4 only
  if (btypes=="mtnd4") biochem = gmm_cluster(ctrls,pat,cluster_chans=MTND4)
  # CIII only
  if (btypes=="cIII") biochem = gmm_cluster(ctrls,pat,cluster_chans=CIII)
  # MTCYB only
  if (btypes=="mtcyb") biochem = gmm_cluster(ctrls,pat,cluster_chans=CIII)  
  # UQCRC2 only
  if (btypes=="uqcrc2") biochem = gmm_cluster(ctrls,pat,cluster_chans=CIII)
  # CIV only
  if (btypes=="cIV") biochem = gmm_cluster(ctrls,pat,cluster_chans=CIV)
  # MTCO1 only
  if (btypes=="mtco1") biochem = gmm_cluster(ctrls,pat,cluster_chans=MTCO1)
  # COX4 only
  if (btypes=="cox4") biochem = gmm_cluster(ctrls,pat,cluster_chans=COX4)
  # CV only
  if (btypes=="cV") biochem = gmm_cluster(ctrls,pat,cluster_chans=CV)
  # OSCP only
  if (btypes=="oscp") biochem = gmm_cluster(ctrls,pat,cluster_chans=OSCP)
  # ATP8 only
  #if (btypes=="atp8") biochem = gmm_cluster(ctrls,pat,cluster_chans=ATP8)
  # ATP5B only
  if (btypes=="atp5b") biochem = gmm_cluster(ctrls,pat,cluster_chans=ATP5B)
  
  # SIGNALLING only
  # OXPHOS & SIGNALLING
  #if (btypes=="signal") biochem = multiD_cluster(pat, Nclust = 1:2, cluster_chans=justprots[!justprots%in%oxp])
  # OXPHOS & SIGNALLING & MORPHOLOGY
  #if (btypes=="oxsig") biochem = multiD_cluster(pat, Nclust = 3, cluster_chans=justprots)
  #if (btypes=="oxsigmorph") biochem = multiD_cluster(pat, Nclust = 3, cluster_chans=c(justprots,morph))
  # Htra2 and HKII
  #if (btypes=="sig1") biochem = multiD_cluster(pat, Nclust = 1:2, cluster_chans=c(Sig1))
  # Linear trend sig with small increase at end
  #if (btypes=="sig2") biochem = multiD_cluster(pat, Nclust = 1:2, cluster_chans=c(Sig2))
  # Unfolded protein response
  #if (btypes=="upr") biochem = multiD_cluster(pat, Nclust = 1:2, cluster_chans=c(UPR))
  # Biogenesis signalling
  #if (btypes=="biogen") biochem = multiD_cluster(pat, Nclust = 1:2, cluster_chans=c(Biogen))
  # Mito protein synth
  #if (btypes=="mitopro") biochem = multiD_cluster(pat, Nclust = 1:2, cluster_chans=c(Mitopro))
  
  # Note that length of pat should be a multiple of the length of biochem:
  length(pat$Value)/length(biochem)
  # There are 58 channels in pat (including LOG versions)
  
  # Add biochem to pat
  # Generate a dummy "channel"
  dummy = pat[pat$Chan=="Area",]
  dummy$Chan = "RCstatus"
  dummy$Channel = "RCstatus"
  dummy$Value = biochem
  dummy$Corrected = biochem
  
  # Patch onto end of pat
  pat = rbind(pat,dummy)
  patlist[[patid]] = pat
}

# Check data are sensible
head(patlist[["P01"]])

pdatnew = do.call("rbind",patlist)
head(pdatnew) 

###Generate plots of clustering
pdf(paste(btypes,"GMM_biochem.pdf",sep=""),width=11.69,height=8.27)
pdat_wide = reshape(pdatnew[,c("Corrected","Chan","ID","SectionID")], idvar = c("ID","SectionID"), 
                    timevar = "Chan", direction="wide")
ctrls_wide = reshape(ctrls[,c("Corrected","Chan","ID","SectionID")], idvar = c("ID","SectionID"),
                    timevar = "Chan", direction="wide")
colnames(pdat_wide) = gsub("Corrected.","",colnames(pdat_wide))
colnames(ctrls_wide) = gsub("Corrected.","",colnames(ctrls_wide))

pats = sort(unique(pdat_wide$SectionID))
for(pat in pats){
  mitochan = "VDAC1"
  df = subset(pdat_wide, SectionID==pat)
  allchans = gsub("LOG_","",colnames(dat)[grepl("LOG",colnames(dat))])
  chans = allchans
  chanvec = chans[chans%in%colnames(df)]
  
  ##Calculate percentages
  N = length(df$ID)
  Pos = sum(df$RCstatus=="Normal",na.rm=TRUE)/N
  nk = sum(df$RCstatus=="Unknown",na.rm=TRUE)/N
  Def = sum(df$RCstatus=="Deficient",na.rm=TRUE)/N
  mlab = paste0("Pos: ",signif(100*Pos,3),"% Def: ", signif(100*Def,3), "%")

  cols = c(rgb(1,0,0,0.2), rgb(0.6,0.2,0.8,0.2), rgb(0,0,1,0.2)) 
  names(cols) = c("Deficient","Unknown","Normal")
  
  # multiD_cluster does not define fibres as "Normal", "Deficient", "Unknown", since clusters don't have that oxphos defect interpretation
  # So, for output from multiD_cluster, we should use a different colouring scheme
  if(typeof(df$RCstatus)=="double"){
    cols = rainbow(n=max(df$RCstatus),alpha=0.2)

    # Plot label which doesn't refer to oxphos status
    ntab = 100*table(df$RCstatus)/N
    mlab = paste(paste0(paste0("C",names(ntab),": "),signif(ntab,3),"%"),collapse=" ") 
  }

  df$col = cols[df$RCstatus] 
  ##Plotting
  op = par(mfrow=c(3,4),mai=c(0.55,0.6,0.2,0.1))  
  for(chan in chans[chans!=mitochan]){
   xpat = log(as.numeric(df[[mitochan]]))
   ypat = log(as.numeric(df[[chan]]))
   xctrl = log(as.numeric(ctrls_wide[[mitochan]]))
   yctrl = log(as.numeric(ctrls_wide[[chan]]))

   # https://stackoverflow.com/questions/16225530/contours-of-percentiles-on-level-plot
   dens = kde2d(xctrl, yctrl, n=200); ## estimate the z counts

   prob = c(0.95, 0.5)
   dx = diff(dens$x[1:2])
   dy = diff(dens$y[1:2])
   sz = sort(dens$z)
   c1 = cumsum(sz) * dx * dy
   levs = sapply(prob, function(x) {
    approx(c1, sz, xout = 1 - x)$y
   })

   if(sum(is.na(df[[chan]]))==length(df[[chan]])){
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5,chan)
   }else{
    plot(xpat,ypat, pch=19, col=df$col, xlab=paste("LOG",mitochan,sep="_"), ylab = paste("LOG",chan,sep="_"),
	  main = mlab,cex.lab=2,cex.main=1.0,cex.axis=1.75,xlim=range(c(xpat,xctrl)),ylim=range(c(ypat,yctrl), na.rm=TRUE, finite=TRUE))    
    abline(v=1.34, col="black", lty=4, lwd=1) # 1.34 represents the bootstrapped cut of for predicting RRF
    contour(dens, levels=levs, labels=prob, add=T,lwd=2)
    }
  }
  mtext(pat, outer=TRUE,  cex=1, line=-2.0,at=0.52,side=3) ##problem with alpha in line above
  par(op)
}
dev.off()

##Write table with additional column to store data, if rerun clustering rewrite table 
write.csv(pdat_wide, file=paste(btypes,"GMM_biochem.csv",sep=""),row.names=FALSE)
}


##investigating area in sig1 clustering
Data<-read.csv("OS_data.csv", header=T)

Data$Chan = Data$Channel
Data$SectionID = pdat$patient_id
ids = sort(unique(Data$patient_id))
ctrl_ids = c("C01","C02","C03","C04")
patids = ids[!ids%in%ctrl_ids]

#Function SEM
sem = function(x){ 
  sd(x)/sqrt(length(x))
  }

pdf("stripchart_allchan_CIdef_CIIIup.pdf",width=6,height=8)
pats = sort(unique(Data$patient_id))
for(patid in patids){
  df = subset(Data, patient_id==patid)
  Num = length(df$ID)
  op = par(mfrow=c(2,3))
  for(m in gsub("-",".",nolog)){
   my_mean = mean(df[[m]])
   my_sem = sem(df[[m]])
   stripchart(df[[m]] ~ df$RC.status,vertical=TRUE,method="jitter",jitter=0.1,pch=16,col=rgb(1,0,0,0.1),cex=1.25,
   main = paste(patid, "N =", Num), xlab = "Group", ylab = m,cex.lab=1.35)
   boxplot(df[[m]] ~ df$RC.status, data = df,add=TRUE,outline=FALSE,notch=TRUE,col=NULL,boxwex=0.3,lwd=2)
  }
  par(op)
}
dev.off()

stripchart(Area ~ RCstatus, data = df,
            frame = FALSE, vertical = TRUE,
            method = "jitter", pch = c(19, 19),
            col = c(rgb(1, 0, 0, 0.2), rgb(0, 0, 1, 0.2)),
            main = paste("Fibre area", patid, "N=", Num), xlab = "Cluster", ylab = "Fibre Area")