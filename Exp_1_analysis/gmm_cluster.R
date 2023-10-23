library(mclust)
library(Hotelling)

getvec = function(dat){
 pc = princomp(dat)
 vec = pc$loadings[,1]
 return(list(vec=vec,cen=pc$center))
}

gmm_cluster = function(ctrls,sec,pcutoff = 0.6,Nsamps = 150000, cluster_chans = c("MTCO1","VDAC1","NDUFB8"), dist_cutoff = 0.5, cluster_cutoff = 0.95){
    sec = sec[sec$Channel%in%unique(ctrls$Channel),]
    ctrls$Corrected = as.numeric(ctrls$Corrected)
    sec$Corrected = as.numeric(sec$Corrected)

    # GMM clustering
    ctrls_wide = reshape(ctrls[,c("Corrected","Chan","ID","SectionID")], #picks data from all row and selected channels
        idvar = c("ID","SectionID"), #define identifiers of row, ID for fibres and SectionID for sections/regions/observer
        timevar = "Chan", direction="wide")
    colnames(ctrls_wide) = gsub("Corrected.","",colnames(ctrls_wide))

    sec_wide = reshape(sec[,c("Corrected","Chan","ID","SectionID")], #picks data from all row and selected channels
        idvar = c("ID","SectionID"), #define identifiers of row, ID for fibres and SectionID for sections/regions/observer
        timevar = "Chan", direction="wide")
    colnames(sec_wide) = gsub("Corrected.","",colnames(sec_wide))

    
    sec_w = log(sec_wide[,cluster_chans])
    ctrls_w = log(ctrls_wide[,cluster_chans])

    keepsec = !is.na(rowSums(sec_w))
    keepctrls = !is.na(rowSums(ctrls_w))

    sec_w = sec_w[keepsec,]
    ctrls_w = ctrls_w[keepctrls,]

    ctrls_wide = ctrls_wide[keepctrls,]
    sec_wide = sec_wide[keepsec,]

    #mb = Mclust(sec_w,2,modelNames=c("VVV"))
    mb = Mclust(sec_w,2,modelNames=NULL)
    
    #if(max(mb$BIC[1,],na.rm=TRUE)>cluster_cutoff*max(mb$BIC[2,],na.rm=TRUE)){ # If two clusters are not much better than #one cluster, just use one:
    # mb = Mclust(sec_w,1)
    #}

    print(mb$modelName)
    samplerows = function(df,N){
         return(df[sample(nrow(df),N,replace=TRUE),])
    }

      diffdfs = function(df1,df2){
         return(as.numeric(sqrt(rowSums((df1-df2)^2))))
    }
      patvals = sec_w[mb$z[,1]>pcutoff,]
      
     if(mb$G==2){
        # If first principle component of c1 is more similar to ctrl than 1st PC of c2 is to ctrl
        c1 = sec_w[mb$z[,1]>pcutoff,]
        c2 = sec_w[mb$z[,1]<(1-pcutoff),]

        diffslopes = function(ctrls_w,c1,c2){
         ct_pc = getvec(ctrls_w)
         c1_pc = getvec(c1)
         c2_pc = getvec(c2)

         c1dist = dist(rbind(ct_pc$vec,c1_pc$vec))[1]
         c2dist = dist(rbind(ct_pc$vec,c2_pc$vec))[1]

         return(c2dist<c1dist)
       }
        ds = replicate(1000,diffslopes(ctrls_w[sample(nrow(ctrls_w), nrow(ctrls_w), replace=TRUE),],c1[sample(nrow(c1), nrow(c1), replace=TRUE),],c2[sample(nrow(c2), nrow(c2), replace=TRUE),]))
        print(mean(ds))

        cert = 0.7 # Should be between 0.5 and 1.0
	  if((mean(ds)>=cert)|(mean(ds)<=(1-cert))){  # IF we are certain about the label order based on slope
           print("Ordering labels based on 1st PC")
           if(mean(ds)>0.5){probs=mb$z[,2]}else{probs=mb$z[,1]}
        }else{
          print("Ordering labels based on Euclidean distance")
          ctmed = apply(ctrls_w,2,median)
          c1med = apply(c1,2,median)
          c2med = apply(c2,2,median)
          d1 = as.numeric(dist(rbind(c1med,ctmed)))
          d2 = as.numeric(dist(rbind(c2med,ctmed)))
          if(d2<d1){probs=mb$z[,2]}else{probs=mb$z[,1]}
        }
       }else{ 
         ctsA = samplerows(ctrls_w,Nsamps)
         ctsB = samplerows(ctrls_w,Nsamps)
         z1 = samplerows(sec_w[mb$z[,1]>pcutoff,],Nsamps)
         ctrlctrl = diffdfs(ctsA,ctsB)
         ctrlpat = diffdfs(ctsA,z1)
         plot(density(ctrlctrl),main="",xlab="Euclidean distance")
         points(density(ctrlpat),col="red",type="l")
         legend("topright",c("Between control fibres","Between patient and control fibres"),col=c("black","red"),lwd=1)
         
         print(mean(ctrlpat))
         if(mean(ctrlpat)>0.5){probs = 1.0-mb$z[,1]}else{probs = mb$z[,1]}

         res = hotelling.test(ctrls_w,sec_w[mb$z[,1]>pcutoff,],perm=TRUE)
         print(res)
    }

     sec_wide$clust = 0
     sec_wide$clust[probs>pcutoff] = 1
     sec_wide$clust[probs<(1-pcutoff)] = 2
     sec_wide$probs = probs

    sec_wide$biochem = "Unknown"
    sec_wide$biochem[sec_wide$clust==1] = "Normal"
    sec_wide$biochem[sec_wide$clust==2] = "Deficient"
    sec_wide$biochem[sec_wide$clust==3] = "Other"

    if(mb$G==1) {
    sec_wide$biochem = "Unknown"
    sec_wide$clust = 0
    }

    bc = sec_wide$biochem
    names(bc) = sec_wide$ID

    return(bc)
}
