#' This function generates results for original simulation study and sensitivity analysis (SNR=2.93, 1.90 and 1.20) performed in Alam and Staicu (20xx)
#' 
#' @param iter is a scalar that represents the iteration number
#' @param n sample size
#' @param nTEST test data sample size
#' @param base_seed a scalar using that a seed number is set for the simulation
#' @export
sim_iter_slfda<-function(iter,n,nTEST,base_seed=100){
  if(nTEST>n)
    stop("Test data size cannot exceed the sample size")
  
  ####### True parameters settings ######
  #SFbasis<-list("sfourier1"=function(s){
  #  fourier(x=s,nbasis=5,period = 1)[,2]
  #},"sfourier2"=function(s){
  #  fourier(x=s,nbasis=5,period = 1)[,3]
  #},"sfourier3"=function(s){
  #  fourier(x=s,nbasis=5,period = 1)[,4]
  #})
  
  
  #TBasis<-list(list("tb11" = function(T){fourier(x=T,nbasis=3,period=1)[,2]},
  #                  "tb12" = function(T){fourier(x=T,nbasis=3,period=1)[,3]}),
  #             list("tb21" = function(T){fourier(x=T,nbasis=5,period=1)[,4]},
  #                  "tb22" = function(T){fourier(x=T,nbasis=5,period=1)[,5]}),
  #             list("tb31" = function(T){fourier(x=T,nbasis=7,period=1)[,6]},
  #                  "tb32" = function(T){fourier(x=T,nbasis=7,period=1)[,7]})
  #)
  
  # Bi-variate basis
  STBasis<-list("stB1"=function(s,t){
    sapply(s, function(u){
      (sqrt(30/37)*(u*t+u^2+t^2))
    })
  },
  "stB2"=function(s,t){
    sapply(s, function(u){
      sqrt(120)*((u*(t^2))-((u^2)*t))
    })
  })
  
  # Marginal basis
  
  SBasis<-list("sB1"=function(s){
    fourier(x=s,nbasis=5,period = 1)[,2]
  },
  "sB2"=function(s){
    fourier(x=s,nbasis=5,period = 1)[,3]
  }
  )
  
  
  
  # To fixed the points on S
  n_s<-51
  ss<-seq(0,1,length.out = n_s)
  
  # Eigen values, eta_{kl}
  tdcfV<-c(0.80,0.60)#list(c(0.6,0.4),c(0.5,0.3),c(0.25,0.20))
  
  
  # Variance parameters
  sig2_S<-c(0.25,0.20)#,0.10)
  sig2_E<-0.15
  
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
  gdata<-SNFDataBB(argS = ss,TimePoint = TijA,STbasis = STBasis,Sbasis = SBasis,
                   Eta = tdcfV,Sigma2K = sig2_S,Sigma2 = sig2_E,
                   muF = meanPF,sclF = sFUN,alpF = alFUN)
  
  trainY<-lapply(seq_len(length(gdata$Y)),function(w){gdata$Y[[w]][1:mi[w],]})
  
  NtestY<-lapply(NTest,function(w){gdata$Y[[w]][mi[w]+1,]})
  
  FtestY<-lapply(FTest,function(w){
    nr<-nrow(gdata$Y[[w]])
    gdata$Y[[w]][(nr-length(Tg)+1):nr,]
  })
  
  rm(gdata)
  
  a<-gc(); rm(a)
  
  ## slfda1
  slfda<-skewedFDA(funDATA=trainY,funARG=ss,obsTIME=Tij,ETGrid=tp,DOP=1,KernelF=depan,CV=TRUE,Hgrid=seq(0.02,0.2,0.02),CVThresh = 0.05,PenaltyF=Qpenalty,plfGT=Tg,
                   ES2knots=c(15,10),ES2bs=c("ps","ps"),ES2m=c(2,2),ES2Esp="REML",
                   LPknots=c(15,10),LPbs=c("ps","ps"),LPm=c(2,2),
                   Cov2nbasis=c(15,10),PVE=c(0.95,0.90),PSid=c(NTest,FTest),PredGS=NULL,PredGT=c(NTij,FTij),
                   Prediction=TRUE,parCOMP=FALSE)
  
  # Next Time Point IPE
  nIPE<-mean(sapply(seq_len(length(NTest)),function(i){
    mean((slfda$PredFD[[i]]-NtestY[[i]])^2)
  }))
  
  # Full process IPE
  
  fIPE<-mean(sapply(seq_len(length(FTest)), function(i){
    mean(apply((slfda$PredFD[[length(NTest)+i]]-FtestY[[i]])^2,2,mean))
  }))
  
  Mean<-do.call(cbind,split(slfda$EstParam$X1,slfda$EstParam$Space))
  Scale<-do.call(cbind,split(slfda$EstParam$X2,slfda$EstParam$Space))
  Alpha<-do.call(cbind,split(slfda$EstParam$X3,slfda$EstParam$Space))[1,]
  h<-slfda$h
  
  # Quantile estimation by slfda1
  qprb<-c(0.05,0.1,0.5,0.9,0.95)
  slfda1q<-cbind(do.call(rbind,quantile_slfda(fitOBJ = slfda,Time=Tg,QLevel = qprb)),rep(qprb,times=length(Tg)),iter)
  
  
  # Out sample validation data
  OTij<-lapply(1:nTEST, function(m){sort(sample(TGrid,sample(5:9,1)))})
  OSmi<-sapply(OTij,length)
  OSPTij<-lapply(1:nTEST, function(m){c(OTij[[m]][OSmi[m]],Tg)})
  OSTij<-lapply(1:nTEST, function(m){c(OTij[[m]],Tg)})
  OSOTij<-lapply(1:nTEST, function(m){OTij[[m]][-OSmi[m]]})
  
  
  vdata<-SNFDataBB(argS = ss,TimePoint = OSTij,STbasis = STBasis,Sbasis = SBasis,
                   Eta = tdcfV,Sigma2K = sig2_S,Sigma2 = sig2_E,
                   muF = meanPF,sclF = sFUN,alpF = alFUN)$Y
  
  Vdata<-lapply(1:nTEST, function(m){vdata[[m]][1:(OSmi[m]-1),]})
  
  TSVdata<-lapply(1:nTEST, function(m){vdata[[m]][-c(1:(OSmi[m]-1)),]})
  
  rm(vdata)
  a<-gc(); rm(a)
  
  pred_slfda<-predict_slfda(fitOBJ = slfda,nfunDATA = Vdata,nfunARG = ss,
                            nobsTIME = OSOTij,identical.ARG = TRUE,PredGT = OSPTij)
  
  ospf1<-mean(sapply(1:nTEST, function(i){
    mean((pred_slfda[[i]][1,]-TSVdata[[i]][1,])^2)
  }))
  
  ospf1f<-mean(sapply(1:nTEST, function(i){
    mean(apply((pred_slfda[[i]][-1,]-TSVdata[[i]][-1,])^2,1,mean))
  }))
  
  
  rm(slfda)
  
  a<-gc(); rm(a)
  # Slfda0
  
  slfda0<-skewedFDA(funDATA=trainY,funARG=ss,obsTIME=Tij,ETGrid=tp,DOP=0,KernelF=depan,CV=TRUE,Hgrid=seq(0.02,0.2,0.02),CVThresh = 0.05,PenaltyF=Qpenalty,plfGT=Tg,
                    ES2knots=c(15,10),ES2bs=c("ps","ps"),ES2m=c(2,2),ES2Esp="REML",
                    LPknots=c(15,10),LPbs=c("ps","ps"),LPm=c(2,2),
                    Cov2nbasis=c(15,10),PVE=c(0.95,0.90),PSid=c(NTest,FTest),PredGS=NULL,PredGT=c(NTij,FTij),
                    Prediction=TRUE,parCOMP=FALSE)
  
  # Next Time Point IPE
  nIPE0<-mean(sapply(seq_len(length(NTest)),function(i){
    mean((slfda0$PredFD[[i]]-NtestY[[i]])^2)
  }))
  
  # Full process IPE
  
  fIPE0<-mean(sapply(seq_len(length(FTest)), function(i){
    mean(apply((slfda0$PredFD[[length(NTest)+i]]-FtestY[[i]])^2,2,mean))
  }))
  
  Mean0<-do.call(cbind,split(slfda0$EstParam$X1,slfda0$EstParam$Space))
  Scale0<-do.call(cbind,split(slfda0$EstParam$X2,slfda0$EstParam$Space)) 
  Alpha0<-do.call(cbind,split(slfda0$EstParam$X3,slfda0$EstParam$Space))[1,]
  h0<-slfda0$h
  
  # Quantile  by slfda0
  slfda0q<-cbind(do.call(rbind,quantile_slfda(fitOBJ = slfda0,Time=Tg,QLevel = qprb)),rep(qprb,times=length(Tg)),iter)
  
  
  pred_slfda0<-predict_slfda(fitOBJ = slfda0,nfunDATA = Vdata,nfunARG = ss,
                             nobsTIME = OSOTij,identical.ARG = TRUE,PredGT = OSPTij)
  
  ospf0<-mean(sapply(1:nTEST, function(i){
    mean((pred_slfda0[[i]][1,]-TSVdata[[i]][1,])^2)
  }))
  
  ospf0f<-mean(sapply(1:nTEST, function(i){
    mean(apply((pred_slfda0[[i]][-1,]-TSVdata[[i]][-1,])^2,1,mean))
  }))
  
  
  rm(slfda0)
  a<-gc(); rm(a)
  
  # Quantile estimation by the COBS
  qprb<-c(0.05,0.1,0.5,0.9,0.95)
  
  # Quantile estimation by COBS
  TrainY<-do.call(rbind,trainY)
  OTij<-do.call(c,Tij)
  cobsQ<-cbind(do.call(rbind,lapply(qprb,function(u){
    NParQ<-sapply(1:n_s,function(v){
      fitC<-cobs(x=OTij,y = TrainY[,v],tau=u,print.mesg = FALSE,print.warn = FALSE)
      predict(fitC,z=Tg)[,2]
    })
    
    Qy<-as.numeric(NParQ)
    Space<-rep(ss,each=length(Tg))
    Time<-rep(Tg,times=length(ss))
    qdat<-data.frame("Qy"=Qy,"Space"=Space,"Time"=Time)
    gfit<-bam(Qy~te(Space,Time,k=c(15,10),bs=c("ps","ps"),m=c(2,2)),data = qdat,method = "REML")
    
    SNParQ<-matrix(gfit$fitted,ncol = length(ss),nrow=length(Tg),byrow = FALSE)
    
    #SNParQ<-sapply(1:nrow(NParQ), function(v){
    #  gam(NParQ[v,]~s(ss,bs="ps",m=2,k=25))$fitted
    #})
    cbind(SNParQ,u)
  })),iter)
  #####
  
  ## LFDA 
  fitLFDA<- fpcaLFDA(Y = do.call(rbind,trainY), subject.index = do.call(c,lapply(1:n, function(i){rep(i,mi[i])})),
                     visit.index = do.call(c,lapply(1:n, function(i){1:mi[i]})), obsT = do.call(c,Tij),
                     funcArg = ss,numTEvalPoints = 51,fbps.knots = c(15,10), fbps.p = 3, fbps.m = 2,
                     LongiModel.method='fpca.sc',
                     mFPCA.pve = 0.95, mFPCA.knots = 15, mFPCA.p = 3, mFPCA.m = 2, 
                     sFPCA.pve = 0.90, sFPCA.nbasis = 10)   
  
  rm(trainY)
  
  a<-gc(); rm(a)
  ## Bivariate Mean by LFDA
  TLim<-range(do.call(c,Tij))
  LTg<-seq(TLim[1],TLim[2],length.out=51)
  meanLFDA<-matrix(predict(object = fitLFDA$meanFOBJ, 
                           newdata = data.frame(x = rep(Tg,length(ss)), 
                                                z = rep(ss, each = length(Tg))))$fitted.values, 
                   nrow = length(Tg))
  
  #qprb<-c(0.05,0.1,0.5,0.9,0.95)
  
  # Quantile estimation by LFDA
  K<-length(fitLFDA$mFPCA.evalues)
  xivar<-sapply(1:K, function(w){apply(sapply(fitLFDA$sFPCA.xiHat.bySubj,function(v){v[,w]}),1,var)})
  qntLFDA<-cbind(do.call(rbind,lapply(c(0.05,0.1,0.5,0.9,0.95), function(w){
    qntw<-qnorm(w,mean=0,sd=sqrt(xivar))%*%t(fitLFDA$mFPCA.efunctions)
    TgQ<-apply(qntw,2,function(u){
      approx(x=fitLFDA$visitTime,y=u,xout = Tg,rule = 2)$y
    })
    cbind(meanLFDA+TgQ,w)
  })),iter)
  
  
  
  LPREV<-predict_lfda(lfdaOBJ = fitLFDA,gridT = c(NTij,FTij),gTID=c(NTest,FTest))  
  
  # Next Time Point IPE for LFDA
  LnIPE<-mean(sapply(seq_len(length(NTest)),function(i){
    mean((LPREV$PredFD[[i]]-NtestY[[i]])^2)
  }))
  
  # Full process IPE for LFDA
  
  LfIPE<-mean(sapply(seq_len(length(FTest)), function(i){
    mean(apply((LPREV$PredFD[[length(NTest)+i]]-FtestY[[i]])^2,1,mean))
  }))
  
  rm(LPREV)
  a<-gc(); rm(a)
  
  # Out of sample prediction by LFDA
  
  osLFDA<-predict_lfda(lfdaOBJ = fitLFDA,nfunDATA=Vdata,nfunARG=ss,nobsTIME=OSOTij,
                       identical.ARG=TRUE,gridS=NULL,gridT=OSPTij)
  
  rm(fitLFDA)
  a<-gc(); rm(a)
  
  
  ospfL<-mean(sapply(1:nTEST, function(i){
    mean((osLFDA$PredFD[[i]][1,]-TSVdata[[i]][1,])^2)
  }))
  
  ospfLf<-mean(sapply(1:nTEST, function(i){
    mean(apply((osLFDA$PredFD[[i]][-1,]-TSVdata[[i]][-1,])^2,1,mean))
  }))
  
  # Returning results
  list("Mean"=rbind(cbind(Mean,1),cbind(Mean0,0),cbind(meanLFDA,2)),
       "Scale"=rbind(cbind(Scale,1),cbind(Scale0,0)),
       "Alpha"=rbind(c(Alpha,1),c(Alpha0,0)),
       "IPE"=c(nIPE,fIPE,nIPE0,fIPE0,LnIPE,LfIPE,ospf1,ospf0,ospf1f,ospf0f,ospfL,ospfLf,sd_n),
       "OptimalBW"=h,
       "Iter"=iter,
       "sLFDA1qntl"=slfda1q,
       "sLFDA0qntl"=slfda0q,
       "LFDAqntl"=qntLFDA,
       "COBSqntl"=cobsQ)
  # Prediction results
  # Results for Next time
}

# Simulation for different shape functions

##Required Libraries 
library(tidyverse)
library(dplyr)
library(MASS)
library(parallel)
library(Rmpi) # not needed if not using the cluster computing on distributed memory
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
Nsim<-1000

for(i in 1:3){
  alFUN<-shape_fun[[i]]
  for(l in 1:3){
    n<-sam_sizes[l]
    ##### Cluster Formation ########
    #cls<-parallel::makeCluster(26)
    cls<-makeCluster( (mpi.universe.size()-1) , type='MPI' )
    registerDoSNOW(cls)
    registerDoParallel(cls)
    
    clusterEvalQ(cls,{
      library(fda)
      library(mgcv)
      library(sn)
      library(refund)
      library(mvtnorm)
      library(cobs)
      library(parallel)
      library(doParallel)
      library(doSNOW)
      library(foreach)
    })
    
    clusterExport(cls,c("CVslfda",
                        "fpcaLFDA",
                        "LocLogLikSN",
                        "predict_slfda",
                        "PWPRED",
                        "OutSamPWPRED",
                        "quantile_slfda",
                        "skewedFDA",
                        "SNFData","alFUN","SNFDataBB","predict_lfda"))
    
    a<-0
    j<-1
    repeat{
      b<-(a+25)
      a1<-((j-1)*25)+1
      b1<-(j*25)
      print(c(a,b,a1,b1))
      # Parallel Simulation run
      par_res<-foreach(i=a1:b1,.errorhandling = "remove") %dopar%
        sim_iter_slfda(iter=i,n=n,nTEST=100,base_seed=100)
      
      # Completed number of iterations
      a<-a+length(par_res)
      # update of seed number
      j<-j+1
      
      # Storing results
      plf_dat<-do.call(rbind,lapply(par_res, function(w){
        cbind(rbind(cbind(w$Mean,1),cbind(w$Scale,2),cbind(w$Alpha,3)),
              w$OptimalBW,w$Iter)
      }))
      
      # Directory at NCSU HPC for storing the results
      dirP<-"/share/astaicu/malam3/SkewedFDA/"
      
      write.table(plf_dat,file = paste(dirP,"PLF_",fnam[i],"_MS_",n,".txt",sep=""),
                  append = TRUE,col.names = FALSE,row.names = FALSE)
      
      ipe_dat<-do.call(rbind,lapply(par_res, function(w){
        c(w$IPE,w$OptimalBW,w$Iter)
      }))
      
      write.table(ipe_dat,file = paste(dirP,"IPE_",fnam[i],"_MS_",n,".txt",sep=""),
                  append = TRUE,col.names = FALSE,row.names = FALSE)
      
      slfda1Q<-do.call(rbind,lapply(par_res, function(w){
        w$sLFDA1qntl
      }))
      
      write.table(slfda1Q,file = paste(dirP,"QsLFDA1_",fnam[i],"_MS_",n,".txt",sep=""),
                  append = TRUE,col.names = FALSE,row.names = FALSE)
      
      slfda0Q<-do.call(rbind,lapply(par_res, function(w){
        w$sLFDA0qntl
      }))
      
      write.table(slfda0Q,file = paste(dirP,"QsLFDA0_",fnam[i],"_MS_",n,".txt",sep=""),
                  append = TRUE,col.names = FALSE,row.names = FALSE)
      
      
      lfdaQ<-do.call(rbind,lapply(par_res, function(w){
        w$LFDAqntl
      }))
      
      write.table(lfdaQ,file = paste(dirP,"QLFDA_",fnam[i],"_MS_",n,".txt",sep=""),
                  append = TRUE,col.names = FALSE,row.names = FALSE)
      
      cobsQ<-do.call(rbind,lapply(par_res, function(w){
        w$COBSqntl
      }))
      
      write.table(cobsQ,file = paste(dirP,"QCOBS_",fnam[i],"_MS_",n,".txt",sep=""),
                  append = TRUE,col.names = FALSE,row.names = FALSE)
      
      
      rm(par_res)
      
      # Termination of repeatition
      if(a>=Nsim)
        break
    }
    
    # Stopping the cluster
    stopCluster(cls)
  }
}

mpi.exit()
