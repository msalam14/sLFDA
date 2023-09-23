data(DTI,package="refund") # loading DTI data
DTI<-DTI %>%
  filter(case==1) # selection of MS patients only
ss_cca<-1:93/93 # conversion of tract locations to be in [0,1]
dti_OT<-DTI$visit.time/max(DTI$visit.time) # conversion of visit time to be in [0,1]
dtiYR<-round(DTI$visit.time/365,2) # visit time in year

misINFO<-which(is.na(DTI$cca),arr.ind = TRUE)
misFA<-split(misINFO[,2],misINFO[,1])
misFP<-as.numeric(names(misFA))
misTRCT<-lapply(misFA,as.numeric)
yDTI<-DTI$cca
for(i in seq_len(length(misTRCT))){
  y<-yDTI[misFP[i],-misTRCT[[i]]]
  x<-ss_cca[-misTRCT[[i]]]
  x_new<-ss_cca[misTRCT[[i]]]
  y_new<-spline(x=x,y=y,xout=x_new,method = "natural")$y
  DTI$cca[misFP[i],misTRCT[[i]]]<-y_new
}


# Construction of patient ID starting from 1
dti_ID<-NULL
dti_ID[1]<-1
for(i in 2:nrow(DTI)){
  a<-ifelse((DTI$ID[i]-DTI$ID[(i-1)])!=0,dti_ID[(i-1)]+1,dti_ID[(i-1)])
  dti_ID[i]<-a
}

# Creating a meta data frame for MS patients
ms_meta<-data.frame("ID"=dti_ID,
                    "ObsTime"=dtiYR)

# Data for implementing the asymmetry test
Y<-split.data.frame(x = DTI$cca,f = dti_ID) # splitting the data for every subject
Tij<-split(x=dti_OT,f=dti_ID) # splitting the time for every subject

# Asymmetry test
time=do.call(c,Tij)
nbreaks<-quantile(time[time>0],probs = seq(0,1,0.12))
nbreaks[c(1,length(nbreaks))]<-c(min(time[time>0])-1e-5,1)
nbreaks<-c(0,nbreaks)
set.seed(130)
asym_test_res<-lfd_asym_test(n_breaks=nbreaks,data=do.call(rbind,Y),time=do.call(c,Tij),alpha = 0.05,boot_ss = 100)
str(asym_test_res)
