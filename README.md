# sLFDA Vignette

# Introduction

This R package is built for analyzing longitudinal functional data
analysis accounting the skewness in marginal distribution. The
methodology based on that it was built is available in Alam and Staicu
(20xx). Note that the package is still in developing stage.

# Installation

To install the package, run the following codes

    devtools::install_github("https://github.com/msalam14/sLFDA")

# DTI data analysis

- First call the required packages

``` r
library(sLFDA)
library(refund)
library(fda)
library(mgcv)
library(sn)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(ggplot2)
```

- Data pre-processing

``` r
data(DTI,package="refund") # loading DTI data
DTI<-DTI %>%
  filter(case==1) # selection of MS patients only
ss_cca<-1:93/93 # conversion of tract locations to be in [0,1]
dti_OT<-DTI$visit.time/max(DTI$visit.time) # conversion of visit time to be in [0,1]
dtiYR<-round(DTI$visit.time/365,2) # visit time in year
```

- There are few FA profiles for which, we do not have FA values for all
  $93$ locations. We fit a spline smoothing for those FA profiles to
  impute the missing values.

``` r
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
```

- Original DTI data has patients’ ID starting from $1001$. We convert
  them to $1,2,\ldots$

``` r
# Construction of patient ID starting from 1
dti_ID<-NULL
dti_ID[1]<-1
for(i in 2:nrow(DTI)){
  a<-ifelse((DTI$ID[i]-DTI$ID[(i-1)])!=0,dti_ID[(i-1)]+1,dti_ID[(i-1)])
  dti_ID[i]<-a
}
```

- Grid selection in the time domain

``` r
ntp<-20 # choice of B
tp<-round(matrix(seq(0,1,length.out = ntp),ncol=1),2) # equally spaced grid of length B
Tg<-seq(0,max(DTI$visit.time),5)/max(DTI$visit.time) # grid where population level functions will be estimated
```

- Data preparation to fit *sLFDA* model

``` r
Y<-split.data.frame(x = DTI$cca,f = dti_ID) # splitting the data for every subject
Tij<-split(x=dti_OT,f=dti_ID) # splitting the time for every subject
```

- Model fitting ($sLFDA_1$)
