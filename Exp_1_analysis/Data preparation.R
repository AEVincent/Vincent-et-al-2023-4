library(tidyr)
library(ggplot2)

#d<-read.csv("mitocyto_merged.out", header=T)

pdat = read.delim("mitocyto_merged_results.csv",sep=",",stringsAsFactors=FALSE)


# Generate a vector of IDs
pats = c(sprintf("P%02d",1:12),sprintf("C%02d",1:4))

# Name the elements of the vector using key from separate file.  This links M-numbers with new patient IDs "M0670-16","M0675-16"
names(pats) = c("M0664-16","M0666-16","M0668-16","M0667-16","M0665-16","M0164-12","M0167-16","M0674-16", "M0670-16", "M0675-16", "M1105-10","M0109-19","M1180-16","M1212-16","M1217-16","M1089-16")

pdat = pdat[pdat$Filename %in% names(pats),]

# Drop P08 because patient is missing
#pats = pats[pats!= "P12"]
#pats = pats[pats!="C03"]

pdat = pdat[pdat$Filename %in% names(pats),]

# Generate patient_id from named vector linking M-numbers with patient IDs:
pdat$patient_id = pats[pdat$Filename]

# Same trick to rename channels: set up a named vector and use it to map from ugly channel names to pretty channel names
chans = c("MTCYB", "LOG_MTCYB", "VDAC1", "LOG_VDAC1", "SDHA", "LOG_SDHA", "MTCO1", "LOG_MTCO1", "NDUFB8", "LOG_NDUFB8", "OSCP", "LOG_OSCP", "MTND4", "LOG_MTND4", "COX4", "LOG_COX4", "ATP5B", "LOG_ATP5B",
          "UqCRC2", "LOG_UqCRC2", "Dystrophin", "LOG_Dystrophin", "Mean","Edges","Area","AspectRatio","Perimeter","Circularity","xCoord","yCoord")

names(chans) = c("115In_CYB", "LOG_115In_CYB", "147Sm_VDAC1", "LOG_147Sm_VDAC1", "153Eu_SDHA", "LOG_153Eu_SDHA", "158Gd_MTCO1", "LOG_158Gd_MTCO1", "160Gd_NDUFB8", "LOG_160Gd_NDUFB8", "161Dy_OSCP", "LOG_161Dy_OSCP", "164Dy_ND4", "LOG_164Dy_ND4", "168Er_COX4", "LOG_168Er_COX4", "170Yb_ATP5B", "LOG_170Yb_ATP5B",
                 "174Yb_UqCRC2", "LOG_174Yb_UqCRC2", "176Yb_Dystrophin", "LOG_176Yb_Dystrophin", "Mean","Edges","Area","AspectRatio","Perimeter","Circularity","xCoord","yCoord")

# Discard data where channel name not listed above
pdat = pdat[pdat$Channel %in% names(chans),]

# Then rename as above
pdat$Channel = chans[pdat$Channel]

# Named vectors are really useful, have a look at them:
print(chans)
print(pats)

###reshape data to wide
dat <- spread(pdat, Channel, Value)
head(dat)

###Generate control df
datac = subset(dat, patient_id%in%c("C01","C02","C03","C04"))

