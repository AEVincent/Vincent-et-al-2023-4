###Read data file
#d<-read.csv("merged_pixel_GMM.csv", header=T)
data = read.delim("merged_GMM_RRF.csv",sep=",",stringsAsFactors=FALSE) #Single pixel data

RRFdat=data%>%filter(Manual_RRF==TRUE)
  
###Define RC conditions 
# NOTE: & not && and "CIV.status" not CIVstatus: would be much better if columns were named systematically

#condition1 = data$CIStatus=="Deficient" & data$CIIIStatus=="Normal" & data$"CIV.status"=="Normal"
#condition2 = data$CIStatus=="Normal" & data$CIIIStatus=="Deficient" & data$"CIV.status"=="Normal"
#condition3 = data$CIStatus=="Normal" & data$CIIIStatus=="Normal" & data$"CIV.status"=="Deficient"
#condition4 = data$CIStatus=="Deficient" & data$CIIIStatus=="Deficient" & data$"CIV.status"=="Normal"
#condition5 = data$CIStatus=="Deficient" & data$CIIIStatus=="Normal" & data$"CIV.status"=="Deficient"
#condition6 = data$CIStatus=="Normal" & data$CIIIStatus=="Deficient" & data$"CIV.status"=="Deficient"
#condition7 = data$CIStatus=="Deficient" & data$CIIIStatus=="Deficient" & data$"CIV.status"=="Deficient"

####Logical function for classification of RC complex deficiency

# NOTE: When testing for equality (i.e. asking a question, not making a statement) use "==" not "="
# NOTE: Could replace "condition1 == TRUE" with just "condition1".  Very useful to think about why that is. 
# NOTE: "class" is a reserved function/keyword in R.  Should not attempt to over-ride with same named function.
# NOTE: What is x in this function?  It is not used at all in the body.  There is no link between the function arguments and the body.
# NOTE: What is returned by this function?  It has only side-effects, printing stuff to the screen.

#classification = function(x){
#  if(condition1 == TRUE){ 
#    print("CI")
#  } 
#  else if(condition2 == TRUE){
#    print("CIII")
#  } 
#  else if(condition3 == TRUE){
#    print("CIV")
#  } 
#  else if(condition4 == TRUE){
#    print("CI+CIII")
#  } 
#  else if(condition5 == TRUE){
#    print("CI+CIV")
#  } 
#  else if(condition6 == TRUE){
#    print("CIII+CIV")
#  } 
#  else if(condition7 == TRUE){
#    print("All")
#  }else {
#  print("None")
#  }
#}

# Seems to be an error in column names:
#colnames(data) = gsub("CiStatus","CIStatus",colnames(data))

###All possible combinations of deficiency
# NOTE: The approach in the following few lines is more concise and works
#data["CI"] = data$CIStatus=="Deficient" & data$CIIIStatus=="Normal" & data$CIVStatus=="Normal"
#data["CIII"] = data$CIStatus=="Normal" & data$CIIIStatus=="Deficient" & data$CIVStatus=="Normal"
#data["CIV"] = data$CIStatus=="Normal" & data$CIIIStatus=="Normal" & data$CIVStatus=="Deficient"
#data["CI+CIII"] = data$CIStatus=="Deficient" & data$CIIIStatus=="Deficient" & data$CIVStatus=="Normal"
#data["CI+CIV"] = data$CIStatus=="Deficient" & data$CIIIStatus=="Normal" & data$CIVStatus=="Deficient"
#data["CIII+CIV"] = data$CIStatus=="Normal" & data$CIIIStatus=="Deficient" & data$CIVStatus=="Deficient"
#data["All"] = data$CIStatus=="Deficient" & data$CIIIStatus=="Deficient" & data$CIVStatus=="Deficient"
#data["None"] = data$CIStatus=="Normal" & data$CIIIStatus=="Normal" & data$CIVStatus=="Normal"

##Without CIII
#data["CI"] = data$CIStatus=="Deficient" & data$CIVStatus=="Normal"
#data["CIV"] = data$CIStatus=="Normal" & data$CIVStatus=="Deficient"
#data["CI+CIV"] = data$CIStatus=="Deficient" & data$CIVStatus=="Deficient"
#data["None"] = data$CIStatus=="Normal" & data$CIVStatus=="Normal"


###Just Summary of Each complex without isolated deficiencies
#data["CI"] = data$CIStatus=="Deficient"
#data["CIII"] = data$CIIIStatus=="Deficient"
#data["CIV"] = data$CIVStatus=="Deficient"

###Using 95%PI to remove upregulated fibres from deficient cluster
#data["CI"] = data$CI_Status=="Deficient" & data$CIII_Status=="Normal" & data$CIV_Status=="Normal" & data$NDUFB8_over_exp=="FALSE"
#data["CIII"] = data$CI_Status=="Normal" & data$CIII_Status=="Deficient" & data$CIV_Status=="Normal" & data$MTCO1_over_exp=="FALSE"
#data["CIV"] = data$CI_Status=="Normal" & data$CIII_Status=="Normal" & data$CIV_Status=="Deficient" & data$CYB_over_exp=="FALSE"
#data["CI+CIII"] = data$CI_Status=="Deficient" & data$CIII_Status=="Deficient" & data$CIV_Status=="Normal" & data$NDUFB8_over_exp=="FALSE" & data$CYB_over_exp=="FALSE"
#data["CI+CIV"] = data$CI_Status=="Deficient" & data$CIII_Status=="Normal" & data$CIV_Status=="Deficient" & data$NDUFB8_over_exp=="FALSE" & data$MTCO1_over_exp=="FALSE"
#data["CIII+CIV"] = data$CI_Status=="Normal" & data$CIII_Status=="Deficient" & data$CIV_Status=="Deficient" & data$MTCO1_over_exp=="FALSE" & data$CYB_over_exp=="FALSE"
#data["All"] = data$CI_Status=="Deficient" & data$CIII_Status=="Deficient" & data$CIV_Status=="Deficient" & data$NDUFB8_over_exp=="FALSE" & data$MTCO1_over_exp=="FALSE" & data$CYB_over_exp=="FALSE"
#data["None"] = data$CI_Status=="Normal" & data$CIII_Status=="Normal" & data$CIV_Status=="Normal"

##EXP1 MITOpanel
data["CV"] = data$ATP5B_status=="Deficient" & data$MTCO1_status=="Normal" & data$NDUFB8_status=="Normal" & data$CIII_status=="Normal"
data["CIV"] = data$ATP5B_status=="Normal" & data$MTCO1_status=="Deficient" & data$NDUFB8_status=="Normal" & data$CIII_status=="Normal"
data["CI"] = data$ATP5B_status=="Normal" & data$MTCO1_status=="Normal" & data$NDUFB8_status=="Deficient" & data$CIII_status=="Normal"
data["CIII"] = data$ATP5B_status=="Normal" & data$MTCO1_status=="Normal" & data$NDUFB8_status=="Normal" & data$CIII_status=="Deficient"
data["CI+CIII"] = data$ATP5B_status=="Normal" & data$MTCO1_status=="Normal" & data$NDUFB8_status=="Deficient" & data$CIII_status=="Normal"
data["CI+CIV"] = data$ATP5B_status=="Normal" & data$MTCO1_status=="Deficient" & data$NDUFB8_status=="Deficient" & data$CIII_status=="Deficient"
data["CI+CV"] = data$ATP5B_status=="Deficient" & data$MTCO1_status=="Normal" & data$NDUFB8_status=="Deficient" & data$CIII_status=="Normal"
data["CIII+CIV"] = data$ATP5B_status=="Normal" & data$MTCO1_status=="Deficient" & data$NDUFB8_status=="Deficient" & data$CIII_status=="Deficient"
data["CIII+CV"] = data$ATP5B_status=="Deficient" & data$MTCO1_status=="Normal" & data$NDUFB8_status=="Normal" & data$CIII_status=="Deficient"
data["CIV+CV"] = data$ATP5B_status=="Deficient" & data$MTCO1_status=="Deficient" & data$NDUFB8_status=="Normal" & data$CIII_status=="Normal"
data["CI+CIII+CIV"] = data$ATP5B_status=="Normal" & data$MTCO1_status=="Deficient" & data$NDUFB8_status=="Deficient" & data$CIII_status=="Deficient"
data["CI+CIV+CV"] = data$ATP5B_status=="Deficient" & data$MTCO1_status=="Deficient" & data$NDUFB8_status=="Deficient" & data$CIII_status=="Normal"
data["CIII+CIV+CV"] = data$ATP5B_status=="Deficient" & data$MTCO1_status=="Deficient" & data$NDUFB8_status=="Normal" & data$CIII_status=="Deficient"
data["All"] = data$ATP5B_status=="Deficient" & data$MTCO1_status=="Deficient" & data$NDUFB8_status=="Deficient" & data$CIII_status=="Deficient"
data["None"] = data$ATP5B_status=="Normal" & data$MTCO1_status=="Normal" & data$NDUFB8_status=="Normal" & data$CIII_status=="Normal"



###Run function for all IDs and generate new column
# NOTE: I would probably leave things as they are above, but if you really want that column, it's probably something like this you want:
dt = data
dt$Status = as.character(sapply(1:dim(dt)[1],function(i) names(dt)[which(as.logical(dt[i,]))]))

#Write table for GMM whole cell
write.table(data, file="GMM_RRF_RRFonly.csv",row.names=FALSE,sep=",")

#summ = aggregate(data[,c("CI","CIII","CIV","CI+CIII","CI+CIV","CIII+CIV","All","None")],by=list(data$SectionID),FUN=mean)
#summ = aggregate(data[,c("CI","CIV","CI+CIV","None")],by=list(data$SectionID),FUN=mean)
#summ = aggregate(data[,c("CI","CIII","CIV")],by=list(data$SectionID),FUN=mean)
summ = aggregate(data[,c("CI","CIII","CIV","CV","CI+CIII","CI+CIV","CI+CV","CIII+CIV","CIII+CV","CIV+CV","CI+CIII+CIV","CI+CIV+CV","CIII+CIV+CV","All","None")],by=list(data$SectionID),FUN=mean)

write.table(summ,file="GMM_SummaryOfDeficiencyCount_RRFonly.csv",row.names=FALSE,sep=",",quote=FALSE)




