#' This function generates results for simulation study on computational time performed in Alam and Staicu (20xx)
#' 
#' @param iter is a scalar that represents the iteration number
#' @param n sample size
#' @param nTEST test data sample size
#' @param base_seed a scalar using that a seed number is set for the simulation
#' @param sig2_E variance of the purely random error
#' @export
sim_iter_slfda<-function(iter,n,nTEST,base_seed=100,p_comp=TRUE,ncr=7){
  if(nTEST>n)
    stop("Test data size cannot exceed the sample size")
  
  ####### True parameters settings ######
  SFbasis<-list("sfourier1"=function(s){
    fourier(x=s,nbasis=5,period = 1)[,2]
  },"sfourier2"=function(s){
    fourier(x=s,nbasis=5,period = 1)[,3]
  },"sfourier3"=function(s){
    fourier(x=s,nbasis=5,period = 1)[,4]
  })
  
  
  TBasis<-list(list("tb11" = function(T){fourier(x=T,nbasis=3,period=1)[,2]},
                    "tb12" = function(T){fourier(x=T,nbasis=3,period=1)[,3]}),
               list("tb21" = function(T){fourier(x=T,nbasis=5,period=1)[,4]},
                    "tb22" = function(T){fourier(x=T,nbasis=5,period=1)[,5]}),
               list("tb31" = function(T){fourier(x=T,nbasis=7,period=1)[,6]},
                    "tb32" = function(T){fourier(x=T,nbasis=7,period=1)[,7]})
  )
  
  
  # To fixed the points on S
  n_s<-51
  ss<-seq(0,1,length.out = n_s)
  
  # Eigen values, eta_{kl}
  tdcfV<-list(c(0.6,0.4),c(0.5,0.3),c(0.25,0.20))
  
  
  # Variance parameters
  sig2_S<-c(0.328,0.210,0.046)
  sig2_E<-0.183
  
  # PLF functions
  # True Mean function
  meanPF<-function(s,t){
    2.5+(3.15*s)+(4*t)+(2*s*t)
  }
  
  #Scale Function
  sFUN<-function(s,t){
    require(mvtnorm)
    25*dmvnorm(c(s,t),mean = c(0,0),sigma=matrix(c(2.5,0.75,0.75,3.5),2,2))
  }
  
  # Testing 
  ntp<-20
  tp<-matrix(seq(0,1,length.out = ntp+2)[-c(1,ntp+2)],ncol=1)
  #' Epanechnikov Kernel function
  #' 
  #' @param x point at which kernel will be evaluated
  #' @export
  depan<-function(x){
    ifelse(abs(x)<=1,(3/4)*(1-x^2),0)
  }
  
  ## Prediction time points
  TGrid<-seq(0,1,length.out=1001)
  
  
  # Full process time grid
  Tg<-seq(0,1,length.out=51)
  ###########################################################
  
  # Specifying iteration specific seed
  sd_n<-((n/base_seed)*1000)+iter
  set.seed(sd_n)
  # Selection of m_i
  mi<-sapply(1:n,sample,x=8:12,size = 1)
  # Subject id for test set for full process and next visit
  FTest<-sample(1:n,nTEST)
  NTest<-sample(1:n,nTEST)
  # Observation time points
  TijA<-lapply(seq_len(length(mi)), function(u){
    if(any(NTest==u) & any(FTest==u)){
      c(sort(sample(TGrid,mi[u]+1)),Tg)
    } else if(any(NTest==u)) {
      sort(sample(TGrid,mi[u]+1))
    } else if(any(FTest==u)){
      c(sort(sample(TGrid,mi[u])),Tg)
    } else{
      sort(sample(TGrid,mi[u]))
    }
  })
  
  # Observed Time points
  Tij<-lapply(seq_len(length(TijA)), function(w){TijA[[w]][1:mi[w]]})
  
  # Next Time point
  NTij<-lapply(NTest, function(w){TijA[[w]][mi[w]+1]})
  
  # Full process time grid
  FTij<-lapply(FTest, function(w){Tg})
  
  # Simulated Skewed FD
  pt<-proc.time()
  gdata<-SNFData(argS = ss,TimePoint = Tij,Sbasis = SFbasis,Tbasis = TBasis,
                 Eta = tdcfV,Sigma2K = sig2_S,Sigma2 = sig2_E,
                 muF = meanPF,sclF = sFUN,alpF = alFUN,parCOMP = FALSE)
  pt0<-proc.time()-pt
  
  trainY<-gdata$Y
  
  rm(gdata)
  
  a<-gc(); rm(a)
  
  # Prediction Timepoints
  predTij<-lapply(Tij, function(t){
    runif(1,max(t),1)
  })
  
  # Computational time assessment
  slfda_ct<-NULL
  
  ## slfda1 with CV
  slfda<-skewedFDA(funDATA=trainY,funARG=ss,obsTIME=Tij,ETGrid=tp,DOP=1,KernelF=depan,CV=TRUE,Hgrid=seq(0.02,0.2,0.02),CVThresh = 0.05,PenaltyF=Qpenalty,plfGT=Tg,
                   ES2knots=c(15,10),ES2bs=c("ps","ps"),ES2m=c(2,2),ES2Esp="REML",
                   LPknots=c(15,10),LPbs=c("ps","ps"),LPm=c(2,2),
                   Cov2nbasis=c(15,10),PVE=c(0.95,0.90),PSid=NULL,PredGS=NULL,PredGT=NULL,
                   Prediction=FALSE,parCOMP = p_comp,n_cores=ncr,fast_bps = FALSE,par_seed=(100+iter))
  slfda_ct<-c(slfda_ct,as.numeric(slfda$CompTime))
  
  # Next time points for the 
  pt_slfda1<-proc.time()
  slfda_pred<-predict_slfda(fitOBJ=slfda,PSid=NULL,PredGS=NULL,PredGT=predTij,nfunDATA=trainY,
                            nfunARG=ss,nobsTIME=Tij,identical.ARG=TRUE,
                            CovDep=FALSE,DesignMat=NULL,NewDesignMat=NULL,PredDesignMat=NULL)
  pt_slfda2<-proc.time()  
  slfda_ct<-c(slfda_ct,as.numeric(pt_slfda2)[3]-as.numeric(pt_slfda1)[3])
  
  # slfda0 with CV
  slfda0<-skewedFDA(funDATA=trainY,funARG=ss,obsTIME=Tij,ETGrid=tp,DOP=0,KernelF=depan,CV=TRUE,Hgrid=seq(0.02,0.2,0.02),CVThresh = 0.05,PenaltyF=Qpenalty,plfGT=Tg,
                    ES2knots=c(15,10),ES2bs=c("ps","ps"),ES2m=c(2,2),ES2Esp="REML",
                    LPknots=c(15,10),LPbs=c("ps","ps"),LPm=c(2,2),
                    Cov2nbasis=c(15,10),PVE=c(0.95,0.90),PSid=NULL,PredGS=NULL,PredGT=NULL,
                    Prediction=FALSE,parCOMP = p_comp,n_cores=ncr,fast_bps = FALSE,par_seed=(100+iter))
  slfda_ct<-c(slfda_ct,as.numeric(slfda0$CompTime))
  
  # Next time points for the 
  pt_slfda1<-proc.time()
  slfda_pred<-predict_slfda(fitOBJ=slfda0,PSid=NULL,PredGS=NULL,PredGT=predTij,nfunDATA=trainY,
                            nfunARG=ss,nobsTIME=Tij,identical.ARG=TRUE,
                            CovDep=FALSE,DesignMat=NULL,NewDesignMat=NULL,PredDesignMat=NULL)
  pt_slfda2<-proc.time()  
  slfda_ct<-c(slfda_ct,as.numeric(pt_slfda2)[3]-as.numeric(pt_slfda1)[3])
  
  
  # lfda
  cp_lfda1<-proc.time()
  fitLFDA<- fpcaLFDA(Y = do.call(rbind,trainY), subject.index = do.call(c,lapply(1:n, function(i){rep(i,mi[i])})),
                     visit.index = do.call(c,lapply(1:n, function(i){1:mi[i]})), obsT = do.call(c,Tij),
                     funcArg = ss,numTEvalPoints = 51,fbps.knots = c(15,10), fbps.p = 3, fbps.m = 2,
                     LongiModel.method='fpca.sc',
                     mFPCA.pve = 0.95, mFPCA.knots = 15, mFPCA.p = 3, mFPCA.m = 2, 
                     sFPCA.pve = 0.90, sFPCA.nbasis = 10)
  cp_lfda2<-proc.time()
  slfda_ct<-c(slfda_ct,as.numeric(cp_lfda2)[3]-as.numeric(cp_lfda1)[3])
  
  cpp_lfda1<-proc.time()
  LPREV<-predict_lfda(lfdaOBJ = fitLFDA,gridS=ss,gTID=NULL,gridT = predTij,
                      nfunDATA=trainY,nfunARG=ss,nobsTIME=Tij,identical.ARG=TRUE) 
  cpp_lfda2<-proc.time()
  slfda_ct<-c(slfda_ct,as.numeric(cpp_lfda2)[3]-as.numeric(cpp_lfda1)[3])
  
  # Returning results
  slfda_ct
}

# Simulation for different shape functions

##Required Libraries 
library(tidyverse)
library(dplyr)
library(MASS)
library(parallel)
library(snow)
library(doSNOW)
library(doParallel)
library(foreach)
library(doRNG)
##################
library(fda)
library(mgcv)
library(sn)
library(refund)
library(mvtnorm)
library(cobs)
##############

#Functions needed for Analysis ####
source("utilities_main.R")

#Shape Function
shape_fun<-list()
shape_fun[[1]]<-function(s){
  2e1*(exp(2*s)/(1+exp(2*s)))*sin(6*pi*s/4)
}
shape_fun[[2]]<-function(s){
  5*(exp(2*s)/(1+exp(2*s)))-0.5*sin(2*pi*s)
}

shape_fun[[3]]<-function(s){
  0
}


fnam<-c("PNS","SN","NRM")
sam_sizes<-c(100,300,500)
Nsim<-100

for(i in 1:3){
  alFUN<-shape_fun[[i]]
  for(l in 1:3){
    n<-sam_sizes[l]
    a<-0
    j<-1
    repeat{
      b<-(a+25)
      a1<-((j-1)*25)+1
      b1<-(j*25)
      print(c(a,b,a1,b1))
      # Parallel Simulation run
      par_res<-foreach(i=a1:b1,.errorhandling = "remove",.combine = 'rbind') %do%
        sim_iter_slfda(iter=i,n=n,nTEST=0,base_seed=100,p_comp=TRUE,ncr=8)
      
      # Completed number of iterations
      a<-a+nrow(par_res)
      # update of seed number
      j<-j+1
      
      # directory in Duke HPC cluster 
      dirP<-"/hpc/home/ma521/"
      
      write.table(par_res,file = paste(dirP,"CompTime_",fnam[i],"_",n,".txt",sep=""),
                  append = TRUE,col.names = FALSE,row.names = FALSE)
      
      rm(par_res)
      
      # Termination of repeatition
      if(a>=Nsim)
        break
    }
  }
}

