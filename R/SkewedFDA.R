#' Longitudinal functional data analysis using sLFDA
#' 
#' Analyze longitudinal functional data using skew-normal distribution for point-
#' wise variation according to the methodlogy proposed in Alam and Staicu(20xx).
#'
#' @param funDATA is a list of functional data observed from different subjects
#' @param funARG functional arguments where the underlying function
#' @param obsTIME is a list of time points for subjects in the study
#' @param ETGrid is a vector of regularly spaced time points
#' @param DOP is a scalar represents degree of polynomial for local fitting for mean and scale function
#' @param KernelF is a character represents the name of the standard kernel to be used local linear maximum likelihood estimation
#' @param h the bandwidth for local MLE
#' @param CV is a logical scalar indicates whether step 1 bandwidth will be determined by cross-validation
#' @param Hgrid a numeric vector of bandwidth values for CV
#' @param CVThresh is a scalar represents the minimum increment for likelihood with change in bandwidth
#' @param PenaltyF is the penalty function for shape parameter of the skewnormal distribution
#' @param plfGT is a vector of time points, taking Cartesian products with funARG, yields grid for predicting mean and scale functions
#' @param ES2knots is a  numeric vector of length two represent the numbers of knots for functional and longitudinal domain required for smoothing at step 2 of the estimation
#' @param ES2bs is a vector of two strings represent names of basis functions in S and T
#' @param ES2m order of penalities to be used in step of PLF estimation
#' @param ES2Esp is a string represents the name of method to be used in step 2 for determining tuning parameter
#' @param LPknots is a numeric vector of length 2 represents the knots for tensor product smoothing of estimated W process in order to demean it
#' @param LPbs is a vector of 2 strings represents the names of basis to be used for tensor product smoothing of estimated W
#' @param LPm is a numeric vector of length 2 represents the order of the basis functions specified in LPbs
#' @param Cov2nbasis is a vector of length 2 represents number of basis functions for 2 applications of fpca.sc for deterimining varphi(s) and psi(t)
#' @param PVE is a vector of length 2 contains the percentage of variation explained required for two fpca applications
#' @param Prediction is a logical scalar indicates whether prediction needs to be make or not; prediction can be made for the subjects in the sample
#' @param PSid is a numeric vector of subject IDs for which prediction will be made
#' @param PredGS is a grid in S for prediction
#' @param PredGT is a list of times points where prediction will be made for different subjects specified in PSid
#' @param CovDep is a logical scalar that represents whether mean depends on covariates or not
#' @param DesignMat is a matrix of baseline covariate for subjects in the data
#' @param CovDepStr is a vecotr of length equals to the number of covariates.
#' Each element can be either 1 or 2; 1 represents coefficient depends on s only, 2 means
#' the coefficient depends on both s and t
#' @param PredDesignMat A list of design matrices for every subject for which the prediction would be made
#' @return a list with following objects:
#' \itemize{
#'  \item funARG : a vector represents where the functions were observed
#'  \item obsTIME : a list of observed time points for different subjects
#'  \item EstParam : a data frame with number of columns equals to 5+number of covariates
#'  \item h : bandwidth of local MLE used in the estimation
#'  \item plfOBJ : fitted objects obtained by mgcv::gam at the 2nd stage of estimation
#'  \item LPMean : mean of the latent process used in the computation of marginal variance
#'  \item lpOBJ : fitted object to smooth the latent W process
#'  \item PLFGS : grid points in functional domain where population level functions were estimated
#'  \item PLFGT : grid points in time domain where population level functions were estimated
#'  \item SBasis : estimated basis function obtained from marginal covariance
#'  \item marMU : estimated mean obtained FPCA of marginal variance
#'  \item TBasis : estimated basis for time-dependent coefficients
#'  \item argTB : time points wehre TBasis were estimated
#'  \item xiMU : estimated means for every FPCA via PACE for time-dependent coefficients
#'  \item EtaKL : estimated eigenvalues for FPCA of every time-dependent coefficients
#'  \item Zeta : predicted scores related to time-dependent coefficients
#'  \item Sigma2K : estimated error variance obtained from FPCA of time-dependent coefficients
#'  \item Sigma2 : estimated error variance for measurement error
#'  \item PredFD : a list of predicted data; elements correspond to subjects
#'  \item PredGS : functional grid where prediction is made
#'  \item PredGT : a list of time points where prediction is made; number of
#'  elements is equal to the number of subjects for them prediction is made
#'  \item PredSUB : subjects id for them prediction is made
#' }
#' @example examples/example_estimation_slfda.R
skewedFDA<-function(funDATA,funARG,obsTIME,ETGrid,DOP=1,KernelF,h=0.05,CV=FALSE,Hgrid=NULL,CVThresh=0.05,PenaltyF=Qpenalty,plfGT=NULL,
                    ES2knots=c(10,10),ES2bs=c("ps","ps"),ES2m=c(2,2),ES2Esp="REML",
                    LPknots=c(25,25),LPbs=c("ps","ps"),LPm=c(2,2),
                    Cov2nbasis=c(20,20),PVE=c(0.95,0.95),
                    Prediction=FALSE,PSid=NULL,PredGS=NULL,PredGT=NULL,
                    CovDep=FALSE,DesignMat=NULL,CovDepStr=NULL,PredDesignMat=NULL){
  if(Prediction & is.null(PredGT))
    stop("Provide timepoints for prediction")
  if(Prediction & is.null(PSid) & length(PredGT)!=length(funDATA))
    stop("Provide prediction time points for subjects in the data")
  if(CV & is.null(Hgrid))
    stop("A vector of bandwidth needs to be provided for CV")
  if(!(DOP==0|DOP==1))
    stop("Set either 1 or 0 for DOP")
  if(Prediction & is.null(PSid))
    print("Prediction is done for all subjects in the data")
#  require(mgcv)
#  require(refund)
#  require(sn)
  n<-length(funDATA)
  mi<-sapply(obsTIME,length)
  subID<-rep(1:n,mi)
  Tij<-do.call(c,obsTIME)
  M<-length(Tij)
  p<-length(funARG)
  B<-length(ETGrid)
  data_mat<-do.call(rbind,funDATA)
  ###
  if(CovDep){
    cds<-c(2,CovDepStr,2,1)
  } else{
    cds<-c(2,2,1)
  }

  # Estimation of PLF
  if(CV){
    cvOBJ<-CVslfda(funDATA = funDATA,funARG = funARG,obsTIME = obsTIME,ETGrid = ETGrid,DOP=DOP,
                   KernelF = KernelF,Hgrid = Hgrid,ES2knots = ES2knots,
                   ES2bs = ES2bs,ES2m = ES2m,ES2Esp = ES2Esp,Thresh = CVThresh,
                   CovDep=CovDep,DesignMat=DesignMat,CovDepStr=CovDepStr)
    h<-cvOBJ$OptimalBW
    fitPLF<-cvOBJ$fitOBJ
    rm(cvOBJ)
  } else{
    XMat<-sapply(ETGrid,function(u){u-Tij})
    IdxM<-apply(XMat,2,function(w){which(abs(w)<=h)})
    ## First stage estimates
    locmleres<-lapply(X = 1:p,FUN = function(u){
      vY<-as.numeric(data_mat[,u])
      t(apply(matrix(1:B),1,FUN=function(v){
        p1<-IdxM[[v]]
        Y<-vY[p1]
        Tm<-Tij[p1]
        if(DOP==0){
          Wt<-(1/h)*KernelF(abs(Tm-ETGrid[v])/h)
          if(CovDep){
            cdim<-ifelse(is.matrix(DesignMat),ncol(DesignMat),1)
            if(cdim==1 & !is.matrix(DesignMat)){
              DesignMat<-matrix(DesignMat,ncol=1)
            }
            Xmat<-cbind(rep(1,length(p1)),DesignMat[p1,])
          } else{
            Xmat<-rep(1,length(p1))
          }
          cnp<-as.numeric(sn.mple(x=Xmat,y=Y,w=Wt,penalty = "Qpenalty",opt.method="BFGS")$cp)
          dpp<-as.numeric(cp2dp(cnp,family = "SN"))
          if(CovDep){
            c(cnp[1:(cdim+1)],log(cnp[(cdim+2)]),dpp[(cdim+3)])
          } else{
            c(cnp[1],log(cnp[2]),dpp[3])
          }
        } else{
          if(CovDep){
            cdim<-ifelse(is.matrix(DesignMat),ncol(DesignMat),1)
            if(cdim==1 & !is.matrix(DesignMat)){
              DesignMat<-matrix(DesignMat,ncol=1)
            }
            Xmat<-cbind(rep(1,length(p1)),XMat[p1,v],DesignMat[p1,])
            in_par<-c(mean(Y),runif(1+cdim),log(sd(Y)),runif(1),runif(1))
            as.numeric(optim(par=in_par,fn=LocLogLikSN,t=ETGrid[v],Tij=Tm,Y=Y,meanX=Xmat,scaleX=Xmat,kernel=KernelF,h=h,penalty=PenaltyF)$par)[-c(2,cdim+4)]
          } else{
            Xmat<-cbind(1,XMat[p1,v])
            in_par<-c(mean(Y),runif(1),log(sd(Y)),runif(1),runif(1))
            as.numeric(optim(par=in_par,fn=LocLogLikSN,t=ETGrid[v],Tij=Tm,Y=Y,meanX=Xmat,scaleX=Xmat,kernel=KernelF,h=h,penalty=PenaltyF)$par)[c(1,3,5)]
          }
        }
      }))
    })
    E1datF<-data.frame(sapply(1:length(cds),function(j){
      as.numeric(sapply(locmleres,function(u){u[,j]}))
    }))
    #######
    rm(locmleres)
    ###########
    E1datF$Space<-rep(funARG,each=B)
    E1datF$Time<-rep(ETGrid,times=p)
    fitPLF<-lapply(1:length(cds), function(u){
      frm<-as.formula(paste("X",u,"~",ifelse(cds[u]==2,
                                             "te(Space,Time,k=ES2knots,bs=ES2bs,m=ES2m)",
                                             "s(Space,k=ES2knots[1],bs=ES2bs[1],m=ES2m[1])"),sep=""))
      gam(frm,data = E1datF,method=ES2Esp)
    })
    rm(E1datF)
  }
  #Estimated Mean and Scale Function
  pdata<-data.frame(Space=rep(funARG,each=length(plfGT)),Time=rep(plfGT,times=p))


  ##### Estimation LP process
  LPdat<-data.frame("Space"=rep(funARG,each=length(Tij)),
                    "Time"=rep(Tij,times=p))

  fitPAR<-sapply(seq_len(length(fitPLF)),function(u){
    predict(fitPLF[[u]],newdata=LPdat)
  })

  if(CovDep){
    fullDM<-cbind(1,apply(as.matrix(DesignMat),2,rep,times=p))
  } else{
    fullDM<-1
  }


  lpred<-do.call(cbind,split(rowSums(as.matrix(fitPAR[,1:(length(cds)-2)])*fullDM),LPdat$Space))
  spred<-do.call(cbind,split(exp(fitPAR[,(length(cds)-1)]),LPdat$Space))
  apred<-as.numeric(do.call(cbind,split(fitPAR[,length(cds)],LPdat$Space))[1,])

  rm(fullDM)

  # Estimation of U process
  sdata_mat<-(data_mat-lpred)/spred
  rm(data_mat); rm(lpred); rm(spred)
  estU<-do.call(cbind,lapply(X=1:p, FUN=function(x){
    psn(sdata_mat[,x],dp=cp2dp(c(0,1,as.numeric(dp2cp(c(0,1,apred[x]),family="SN"))[3]),family = "SN"))
  }))

  estU[estU<1e-3]<-1e-3
  estU[estU>(1-1e-3)]<-(1-1e-3)

  # Estimation of W process
  estW<-qnorm(estU)

  rm(estU)

  # Estimation of basis function for functional argument
  # Data preparation
  LPdat$W=as.numeric(estW)
  bam_fit<-bam(W~te(Space,Time,k=LPknots,bs=LPbs,m=LPm),data=LPdat)
  rm(LPdat)
  # Functional principal component analysis
  sc_fit1<-fpca.sc(Y=matrix(bam_fit$residuals,ncol = p,byrow = FALSE),argvals = funARG,nbasis = Cov2nbasis[1],pve = PVE[1],var = TRUE)
  phiH<-sc_fit1$efunctions
  sc1mu<-sc_fit1$mu
  Sigma2<-sc_fit1$sigma2
  rm(sc_fit1)
  psiSC<-((matrix(bam_fit$residuals,ncol = p,byrow = FALSE)-outer(rep(1,M),sc1mu))%*%phiH)*(1/p)
  Khat<-ncol(psiSC)
  TFPCA<-lapply(seq_len(Khat), function(u){
    ydata_k<-data.frame(subID,Tij,psiSC[,u])
    colnames(ydata_k)<-c(".id",".index",".value")
    fpca.sc(ydata = ydata_k,var = TRUE,nbasis = Cov2nbasis[2],pve=PVE[2])
  })

  # Prediction
  PredFD<-NULL
  if(Prediction){
    pmi<-sapply(PredGT,length)
    PTij<-do.call(c,PredGT)
    if(is.null(PSid)){
      PSid<-1:n
    }
    PrID<-rep(PSid,pmi)
    if(is.null(PredGS)){
      PredGS<-funARG
    }
    PredPHI<-sapply(seq_len(Khat),function(w){
      approx(x=funARG,y=phiH[,w],xout = PredGS,rule = 2)$y
    })
    Palp<-predict(fitPLF[[length(cds)]],newdata=data.frame(Space=PredGS))
    Psc1mu<-approx(x=funARG,y=sc1mu,xout = PredGS,rule = 2)$y

    psiKL<-lapply(1:Khat, function(k){
      Lk<-ncol(TFPCA[[k]]$efunctions)
      lapply(1:Lk, function(l){
        approxfun(x=TFPCA[[k]]$argvals,y=TFPCA[[k]]$efunctions[,l],rule = 2)
      })
    })

    xiMU<-lapply(1:Khat, function(k){
      approxfun(x=TFPCA[[k]]$argvals,y=TFPCA[[k]]$mu,rule = 2)
    })


    BVmean<-matrix(predict(bam_fit,newdata=data.frame(Space=rep(PredGS,each=length(PTij)),Time=rep(PTij,times=length(PredGS)))),nrow=length(PTij),byrow=FALSE)
    UstM<-pnorm(Reduce(`+`,lapply(seq_len(Khat), function(k){
      Lk<-length(psiKL[[k]])
      Reduce(`+`,lapply(seq_len(Lk),function(l){
        outer((TFPCA[[k]]$scores[PrID,l]*psiKL[[k]][[l]](PTij)),PredPHI[,k])
      }))+outer(xiMU[[k]](PTij),PredPHI[,k])
    }))+(outer(rep(1,length(PTij)),Psc1mu))+BVmean)

    PSInd<-rep(1:length(PSid),pmi)
    UstM[UstM<1e-3]<-1e-3
    UstM[UstM>(1-1e-3)]<-(1-1e-3)

    if(CovDep){
      if(length(cds)==4){
        PredDesignMat<-do.call(c,PredDesignMat)
      } else{
        PredDesignMat<-do.call(rbind,PredDesignMat)
      }
      fullDM<-cbind(1,apply(as.matrix(PredDesignMat),2,rep,times=length(PredGS)))
    } else{
      fullDM<-1
    }

    predPAR<-sapply(seq_len(length(fitPLF)),function(u){
      predict(fitPLF[[u]],newdata=data.frame("Space"=rep(PredGS,each=length(PTij)),
                                             "Time"=rep(PTij,times=length(PredGS))))
    })

    Lpred<-do.call(cbind,split(rowSums(as.matrix(predPAR[,1:(length(cds)-2)])*fullDM),rep(PredGS,each=length(PTij))))
    Spred<-do.call(cbind,split(exp(predPAR[,(length(cds)-1)]),rep(PredGS,each=length(PTij))))


    PredFD<-split.data.frame((Lpred+(
      Spred*
        sapply(seq_len(ncol(UstM)),function(a){
          qsn(UstM[,a],dp=cp2dp(c(0,1,as.numeric(dp2cp(c(0,1,Palp[a]),family="SN"))[3]),family = "SN"),solver = "RFB")
        })
    )),PSInd)
    rm(UstM)
  }

  # Return Object
  return(list("funARG"=funARG,
              "obsTIME"=obsTIME,
              "EstParam"=data.frame(sapply(seq_len(length(fitPLF)),function(u){
                predict(fitPLF[[u]],newdata=data.frame("Space"=rep(funARG,each=length(plfGT)),
                                                       "Time"=rep(plfGT,times=p)))
              }),"Space"=rep(funARG,each=length(plfGT)),"Time"=rep(plfGT,times=p)),
              "h"=h,
              "plfOBJ"=fitPLF,
              "LPMean"=do.call(cbind,split(predict(bam_fit,newdata=pdata),pdata$Space)),
              "lpOBJ" = bam_fit,
              "PLFGS" = funARG,
              "PLFGT" = plfGT,
              "SBasis"=phiH,
              "marMU" = sc1mu,
              "TBasis"=lapply(TFPCA,function(w){w$efunctions}),
              "argTB"=sapply(TFPCA,function(w){w$argvals}),
              "xiMU" = sapply(TFPCA,function(w){w$mu}),
              "EtaKL" = lapply(TFPCA,function(w){w$evalues}),
              "Zeta" = lapply(TFPCA,function(w){w$scores}),
              "Sigma2K" = sapply(TFPCA,function(w){w$sigma2}),
              "Sigma2" = Sigma2,
              "PredFD"=PredFD,
              "PredGS" = PredGS,
              "PredGT" = PredGT,
              "PredSUB"= PSid))
}


#' Predicts trajectory by the sLFDA method
#' 
#' Performs prediction using the fitted slfda object for subjects no matter whether used in the model fitting or not. However, prediction can be made one type of subjects at a time.
#'
#' @param fitOBJ fitted object obtained through slfda function
#' @param PSid is a numeric vector of subject IDs for which prediction will be made; needed for in-sample prediction only
#' @param PredGS is a grid in S where prediction will be made prediction
#' @param PredGT is a list of times points where prediction will be made; for in-sample prediction the length should be the length of PSid, however,
#' for out-of-sample prediction the length must be equal to the length of the list nfunDATA
#' @param outSAMPLE is a logical scalar indicates whether perdiction is for out-of-sample data or not
#' @param nfunDATA a list containing funcational data for out-of-sample subjects
#' @param nfunARG is a numeric vector contain values where nfunDATA is sampled
#' @param nobsTIME is a list of numeric vectors contains observation time points for subjects in the nfunDATA
#' @param identical.ARG is a logical scalar indicates whether functional data for out-of-sample subjects are have sample at identical sampled points as the functional data used in model fitting
#' @param CovDep is a logical scalar indicates fitOBJ obtained with covariates dependency or not
#' @param DesignMat is a design matrix for the subjects used in the fitOBJ
#' @param NewDesignMat is a design matrices for different subjects that corresponds to the observed time points
#' @param PredDesignMat is a list of design matrices for different subjects that correspond to the prediction time points
#' @return a list with predicted data; elements correspond to subjects
#' @example examples/example_prediction_slfda.R
#' @export
predict_slfda<-function(fitOBJ,PSid=NULL,PredGS=NULL,PredGT=NULL,outSAMPLE=FALSE,
                        nfunDATA=NULL,nfunARG=NULL,nobsTIME=NULL,identical.ARG=TRUE,
                        CovDep=FALSE,DesignMat=NULL,NewDesignMat=NULL,PredDesignMat=NULL){
  funARG<-fitOBJ$PLFGS
  phiH<-fitOBJ$SBasis
  Khat<-ncol(fitOBJ$SBasis)
  sc1mu<-approxfun(x=funARG,y=fitOBJ$marMU,rule = 2)
  psiKL<-lapply(1:Khat, function(k){
    Lk<-ncol(fitOBJ$TBasis[[k]])
    lapply(1:Lk, function(l){
      approxfun(x=fitOBJ$argTB[,k],y=fitOBJ$TBasis[[k]][,l],rule = 2)
    })
  })

  xiMU<-lapply(1:Khat, function(k){
    approxfun(x=fitOBJ$argTB[,k],y=fitOBJ$xiMU[,k],rule = 2)
  })

  if(outSAMPLE){
    if(any(c(is.null(nfunDATA),is.null(nfunARG),is.null(nobsTIME))))
      stop("Out sample prediction requires functional data, argument and observed time points")
    if(!is.null(PSid))
      stop("For out-of-sample prediction, keep PSid as NULL")
    if(any(round(diff(diff(nfunARG)),6)!=0))
      stop("Out sample functions are not observed in fine grid")
    if(identical.ARG & any(round(nfunARG-fitOBJ$PLFGS,6)!=0))
      stop("Check functional argument for identity when identical.ARG is true")
    if(!identical.ARG){
      phiH<-apply(phiH, 2, FUN = function(u){
        approx(x=funARG,y=u,xout = nfunARG,rule=2)$y
      })
    }
    p<-length(nfunARG)
    nID<-do.call(c,lapply(seq_len(length(nobsTIME)), function(w){rep(w,length(nobsTIME[[w]]))}))
    pTij<-do.call(c,nobsTIME)
    data_mat<-do.call(rbind,nfunDATA)
    ##### Estimation LP process
    LPdat<-data.frame("Space"=rep(nfunARG,each=length(pTij)),
                      "Time"=rep(pTij,times=p))

    fitPAR<-sapply(seq_len(length(fitOBJ$plfOBJ)),function(u){
      predict(fitOBJ$plfOBJ[[u]],newdata=LPdat)
    })

    if(CovDep){
      fullDM<-cbind(1,apply(as.matrix(NewDesignMat),2,rep,times=p))
    } else{
      fullDM<-1
    }

    lpred<-do.call(cbind,split(rowSums(as.matrix(fitPAR[,1:(ncol(fitPAR)-2)])*fullDM),LPdat$Space))
    spred<-do.call(cbind,split(exp(fitPAR[,(ncol(fitPAR)-1)]),LPdat$Space))
    if(!identical.ARG){
      apred<-predict(fitOBJ$plfOBJ[[ncol(fitPAR)]],newdata=data.frame(Space=nfunARG))
    } else{
      apred<-as.numeric(do.call(cbind,split(fitPAR[,ncol(fitPAR)],LPdat$Space))[1,])
    }
    # Estimation of U process
    sdata_mat<-(data_mat-lpred)/spred
    estU<-do.call(cbind,lapply(X=1:p, FUN=function(x){
      psn(sdata_mat[,x],dp=cp2dp(c(0,1,as.numeric(dp2cp(c(0,1,apred[x]),family="SN"))[3]),family = "SN"))
    }))

    estU[estU<1e-3]<-1e-3
    estU[estU>(1-1e-3)]<-(1-1e-3)

    # Estimation of W process
    estW<-qnorm(estU)

    rm(estU)

    BVmean<-do.call(cbind,split(predict(fitOBJ$lpOBJ,newdata=LPdat),LPdat$Space))
    dW<-estW-BVmean

    rm(estW); rm(LPdat)
    dSW<-split.data.frame(dW,nID)
    NXi<-do.call(rbind,lapply(seq_len(length(dSW)), function(w){
      if(length(nobsTIME[[w]])==1){
        apply(phiH,2,FUN = function(v){
          sum((dSW[[w]]-sc1mu(nfunARG))*v*diff(c(0,nfunARG)))
        })
      } else{
        t(sapply(seq_len(length(nobsTIME[[w]])), function(j){
          apply(phiH,2,FUN = function(v){
            sum((dSW[[w]][j,]-sc1mu(nfunARG))*v*diff(c(0,nfunARG)))
          })
        }))
      }
    }))
    subGK<-lapply(seq_len(length(nobsTIME)),function(i){
      lapply(1:Khat,function(k){
        Lk<-length(psiKL[[k]])
        Gk<-Reduce(`+`,lapply(1:Lk, function(l){
          fitOBJ$EtaKL[[k]][l]*outer(psiKL[[k]][[l]](nobsTIME[[i]]),psiKL[[k]][[l]](nobsTIME[[i]]))
        }))
        if(fitOBJ$Sigma2K[k]<1e-10){
          solve(Gk+(1e-10*diag(length(nobsTIME[[i]]))))
        } else{
          solve(Gk+(fitOBJ$Sigma2K[k]*diag(length(nobsTIME[[i]]))))
        }
      })
    })
    nZETA<-lapply(seq_len(Khat),function(k){
      knXI<-split(NXi[,k],nID)
      Lk<-length(psiKL[[k]])
      do.call(rbind,lapply(seq_len(length(nobsTIME)),function(i){
        sapply(1:Lk, function(l){
          fitOBJ$EtaKL[[k]][l]*sum(psiKL[[k]][[l]](nobsTIME[[i]])*(subGK[[i]][[k]]%*%matrix(knXI[[i]]-xiMU[[k]](nobsTIME[[i]]),ncol=1)))
        })
      }))
    })
    # Prediction
    if(!is.null(PredGS)){
      if(identical.ARG){
        phiH<-apply(phiH, 2, FUN = function(u){
          approx(x=funARG,y=u,xout = PredGS,rule=2)$y
        })
      } else{
        phiH<-apply(phiH, 2, FUN = function(u){
          approx(x=nfunARG,y=u,xout = PredGS,rule=2)$y
        })
      }
    } else{
      PredGS<-nfunARG
    }

    pgt_ind<-do.call(c,lapply(seq_len(length(PredGT)),function(u){
      rep(u,length(PredGT[[u]]))
    }))
    if(!is.null(PredGT)){
      pTij<-do.call(c,PredGT)
      LPdat<-data.frame("Space"=rep(PredGS,each=length(pTij)),
                        "Time"=rep(pTij,times=p))

      fitPAR<-sapply(seq_len(length(fitOBJ$plfOBJ)),function(u){
        predict(fitOBJ$plfOBJ[[u]],newdata=LPdat)
      })

      if(CovDep){
        if(ncol(fitPAR)==4){
          PredDesignMat<-do.call(c,PredDesignMat)
        } else{
          PredDesignMat<-do.call(rbind,PredDesignMat)
        }
        fullDM<-cbind(1,apply(as.matrix(PredDesignMat),2,rep,times=length(PredGS)))
      } else{
        fullDM<-1
      }

      lpred<-do.call(cbind,split(rowSums(as.matrix(fitPAR[,1:(ncol(fitPAR)-2)])*fullDM),LPdat$Space))
      spred<-do.call(cbind,split(exp(fitPAR[,(ncol(fitPAR)-1)]),LPdat$Space))
      BVmean<-do.call(cbind,split(predict(fitOBJ$lpOBJ,newdata=LPdat),LPdat$Space))
      nID<-do.call(c,lapply(seq_len(length(PredGT)),function(u){
        rep(u,length(PredGT[[u]]))
      }))
    }
    UstM<-pnorm(Reduce(`+`,lapply(seq_len(Khat), function(k){
      Lk<-length(psiKL[[k]])
      Reduce(`+`,lapply(seq_len(Lk),function(l){
        outer((nZETA[[k]][nID,l]*psiKL[[k]][[l]](pTij)),phiH[,k])
      }))+outer(xiMU[[k]](pTij),phiH[,k])
    }))+(outer(rep(1,length(pTij)),sc1mu(nfunARG)))+BVmean)
    UstM[UstM<1e-3]<-1e-3
    UstM[UstM>(1-1e-3)]<-(1-1e-3)
    PredFD<-split.data.frame(lpred+(spred*(
        sapply(seq_len(ncol(UstM)),function(a){
          qsn(UstM[,a],dp=cp2dp(c(0,1,as.numeric(dp2cp(c(0,1,apred[a]),family="SN"))[3]),family = "SN"),solver = "RFB")
        })
    )),nID)
    rm(UstM)
    PredFD
  } else{
    if(any(c(!is.null(nfunDATA),!is.null(nfunARG),!is.null(nobsTIME))))
      stop("Prediction needs to be done for either in-sample or out-sample subjects")
    if(is.null(PSid) & !is.null(PredGT) & length(PredGT)!=nrow(fitOBJ$Zeta[[1]]))
      stop("Provide prediction time points for all subjects in the data used for model fitting")
    if(!is.null(PSid) & length(PredGT)!=length(PSid))
      stop("Provide time points for only the interested subjects")

    if(is.null(PredGT))
      PredGT<-fitOBJ$obsTIME
    pmi<-sapply(PredGT,length)
    PTij<-do.call(c,PredGT)
    n<-length(fitOBJ$obsTIME)
    p<-length(funARG)
    if(is.null(PSid)){
      PSid<-1:n
      if(is.null(PredDesignMat)){
        PredDesignMat<-split.data.frame(as.matrix(DesignMat),rep(PSid,pmi))
      }
    }
    PrID<-rep(PSid,pmi)
    if(is.null(PredGS)){
      PredGS<-funARG
    }
    PredPHI<-sapply(seq_len(Khat),function(w){
      approx(x=funARG,y=phiH[,w],xout = PredGS,rule = 2)$y
    })
    Palp<-predict(fitOBJ$plfOBJ[[length(fitOBJ$plfOBJ)]],newdata=data.frame(Space=PredGS))

    BVmean<-matrix(predict(fitOBJ$lpOBJ,newdata=data.frame(Space=rep(PredGS,each=length(PTij)),Time=rep(PTij,times=length(PredGS)))),nrow=length(PTij),byrow=FALSE)
    UstM<-pnorm(Reduce(`+`,lapply(seq_len(Khat), function(k){
      Lk<-length(psiKL[[k]])
      Reduce(`+`,lapply(seq_len(Lk),function(l){
        outer((fitOBJ$Zeta[[k]][PrID,l]*psiKL[[k]][[l]](PTij)),PredPHI[,k])
      }))+outer(xiMU[[k]](PTij),PredPHI[,k])
    }))+(outer(rep(1,length(PTij)),sc1mu(PredGS)))+BVmean)

    PSInd<-rep(1:length(PSid),pmi)
    UstM[UstM<1e-3]<-1e-3
    UstM[UstM>(1-1e-3)]<-(1-1e-3)

    LPdat<-data.frame("Space"=rep(PredGS,each=length(PTij)),
                      "Time"=rep(PTij,times=p))

    fitPAR<-sapply(seq_len(length(fitOBJ$plfOBJ)),function(u){
      predict(fitOBJ$plfOBJ[[u]],newdata=LPdat)
    })


    if(CovDep){
      if(ncol(fitPAR)==4){
        PredDesignMat<-do.call(c,PredDesignMat)
      } else{
        PredDesignMat<-do.call(rbind,PredDesignMat)
      }
      fullDM<-cbind(1,apply(as.matrix(PredDesignMat),2,rep,times=length(PredGS)))
    } else{
      fullDM<-1
    }

    lpred<-do.call(cbind,split(rowSums(as.matrix(fitPAR[,1:(ncol(fitPAR)-2)])*fullDM),LPdat$Space))
    spred<-do.call(cbind,split(exp(fitPAR[,(ncol(fitPAR)-1)]),LPdat$Space))

    PredFD<-split.data.frame((lpred+(
      spred*
        sapply(seq_len(ncol(UstM)),function(a){
          qsn(UstM[,a],dp=cp2dp(c(0,1,as.numeric(dp2cp(c(0,1,Palp[a]),family="SN"))[3]),family = "SN"),solver = "RFB")
        })
    )),PSInd)
    rm(UstM)
    PredFD
  }
}

#' Estimate quantile trajectory by the sLFDA method
#' 
#' Calculate the quantile trajectory in LFDA using the fitted object obtained by slfda proposed by Alam and Staicu (20XX)
#'
#' @param fitOBJ fitted object obtained through slfda function
#' @param Time is numeric vector where quantile functions will be estimated
#' @param QLevel is a numeric vector represent the level of intended quantile trajectory
#' @param CovDep is a logical scalar indicates whether slfda model was fitted with covariates or not
#' @param NewDesignMat is a design matrix to adjust the mean for covariate effect; only require when CovDep is true
#' @returns a list with predicted quantile trajectories at the given time points
#' @example examples/example_quantile_estimation.R
#' @export
quantile_slfda<-function(fitOBJ,Time,QLevel,
                         CovDep=FALSE,NewDesignMat=NULL){
  funARG<-fitOBJ$PLFGS
  LPdat<-data.frame("Space"=rep(funARG,each=length(Time)),
                    "Time"=rep(Time,times=length(funARG)))

  fitPAR<-sapply(seq_len(length(fitOBJ$plfOBJ)),function(u){
    predict(fitOBJ$plfOBJ[[u]],newdata=LPdat)
  })

  if(CovDep){
    if(ncol(fitPAR)==4){
      NewDesignMat<-do.call(c,NewDesignMat)
    } else{
      NewDesignMat<-do.call(rbind,NewDesignMat)
    }
    fullDM<-cbind(1,apply(as.matrix(NewDesignMat),2,rep,times=length(funARG)))
  } else{
    fullDM<-1
  }

  lpred<-do.call(cbind,split(rowSums(as.matrix(fitPAR[,1:(ncol(fitPAR)-2)])*fullDM),LPdat$Space))
  spred<-do.call(cbind,split(exp(fitPAR[,(ncol(fitPAR)-1)]),LPdat$Space))
  apred<-as.numeric(do.call(cbind,split(fitPAR[,ncol(fitPAR)],LPdat$Space))[1,])

  lapply(1:length(Time), function(t){
    outer(rep(1,length(QLevel)),lpred[t,])+
      outer(rep(1,length(QLevel)),spred[t,])*sapply(
        apred,function(a){qsn(QLevel,dp=cp2dp(c(0,1,dp2cp(c(0,1,a),family="SN")[3]),family="SN"))}
      )
  })
}


#' Skew-normal log-likelihood calculation
#' 
#' Calculate the penalized local skew-normal log-likelihood function where both
#' mean and scale can be approximated by polynomials. This function needs to
#' apply for every functional grid separately when using in slfda.
#'
#' @param par is a vector of parameters; the last element must be the shape parameter, and all the previous are for for local approximations of mean and scale, respectively
#' @param t is a scalar where log-likelihood would be calculated
#' @param Tij is a vector of observed time points stacked for all subjects
#' @param Y is a vector of data corresponding to the time points Tij
#' @param meanX represents design matrix for local approximation of mean and covariate dependency
#' @param scaleX represents design matrix for local approximation of scale
#' @param kernel is a kernel function that calculates weights
#' @param h is a scalar for bandwidth of the given kernel function
#' @param penalty is a function that penalizes the shape parameter
#' @returns a scalar that represents the value of log-likelihood
#' @example examples/example_local_likelihood.R
#' @export
LocLogLikSN<-function(par,t,Tij,Y,meanX,scaleX,kernel,h,penalty=penAZs){
  a<-ncol(meanX)
  b<-ncol(scaleX)
  p<-length(par)
  hij<-abs(Tij-t)
  Wt<-(1/h)*kernel(hij/h)
  ### Calculation of penalty for alpha
  pen<-penalty(par[p])
  ### Parameter for SN likelihood
  bd<-sqrt(2/pi)*(par[p]/sqrt(1+par[p]^2))
  TScale<-exp(as.numeric(scaleX%*%matrix(par[(a+1):(a+b)])))
  TScale<-TScale/sqrt(1-bd^2)
  Loc<-(as.numeric(meanX%*%as.matrix(par[1:a])))-(TScale*bd)
  ### Penalized weighted log-likelihood (negative)
  -(sum(dsn(x=Y,xi=Loc,omega = TScale,alpha = par[p],log = TRUE)*Wt)-pen)
}



#' Cross-validation in two-step estimation of PLF by the sLFDA
#' 
#' Cross-Validation based on likelihood function for bandwidth selection at the
#' first step of estimation proposed for population level function estimation in
#' Alam and Staicu (20xx)
#'
#' @param funDATA is a list of functional data observed from different subjects
#' @param funARG functional arguments where the underlying function
#' @param obsTIME is a list of time points for subjects in the study
#' @param ETGrid is a vector of regularly spaced time points
#' @param DOP is a scalar represents degree of polynomial for local fitting for mean and scale function
#' @param KernelF is a character represents the name of the standard kernel to be used local linear maximum likelihood estimation
#' @param Hgrid a numeric vector of bandwidth values for which log-likelihood will be evaluated
#' @param CVThresh is a scalar represents the minimum increment for likelihood with change in bandwidth
#' @param PenaltyF is the penalty function for shape parameter of the skewnormal distribution
#' @param ES2knots is a  numeric vector of length two represent the numbers of knots for functional and longitudinal domain required for smoothing at step 2 of the estimation
#' @param ES2bs is a vector of two strings represent names of basis functions in S and T
#' @param ES2m order of penalities to be used in step of PLF estimation
#' @param ES2Esp is a string represents the name of method to be used in step 2 for determining tuning parameter
#' @param Thresh is a scalar represents the minimum increment for likelihood with change in bandwidth
#' @param CovDep is a logical scalar that represents whether mean depends on covariates or not
#' @param DesignMat is a matrix of baseline covariate for subjects in the data
#' @param CovDepStr is a vecotr of length equals to the number of covariates.
#' Each element can be either 1 or 2; 1 represents coefficient depends on s only, 2 means
#' the coefficient depends on both s and t
#' @return a list with the following items:
#' \itemize{
#'  \item CVresults : a data frame with two columns (bandwidth and log-likelihood)
#'  \item OptimalBW : chosen optimal bandwidth value
#'  \item fitOBJ : smooth fitted objects for population level functions at the optimal bandwidth
#' }
#' @example examples/exampl_cv_slfda.R
#' @export
CVslfda<-function(funDATA,funARG,obsTIME,ETGrid,DOP=1,KernelF,Hgrid,PenaltyF=Qpenalty,
                  ES2knots=c(10,10),ES2bs=c("ps","ps"),ES2m=c(2,2),ES2Esp="REML",Thresh=0.05,
                  CovDep=FALSE,DesignMat=NULL,CovDepStr=NULL){
  #  require(mgcv)
  #  require(refund)
  #  require(sn)
  if(!(DOP==0|DOP==1))
    stop("Set either 1 or 0 for DOP")
  n<-length(funDATA)
  mi<-sapply(obsTIME,length)
  subID<-rep(1:n,mi)
  Tij<-do.call(c,obsTIME)
  p<-length(funARG)
  B<-length(ETGrid)
  #
  if(CovDep){
    cds<-c(2,CovDepStr,2,1)
  } else{
    cds<-c(2,2,1)
  }

  # Estimation of PLF
  data_mat<-do.call(rbind,funDATA)
  XMat<-sapply(ETGrid,function(u){u-Tij})
  ## First stage estimates
  i<-1
  Nh<-length(Hgrid)
  RC<-1
  Lh<-NULL
  bandw<-NULL
  while(i<=Nh & RC>Thresh){
    h<-Hgrid[i]
    bandw[i]<-h
    IdxM<-apply(XMat,2,function(w){which(abs(w)<=h)})
    locmleres<-lapply(X=1:p,FUN=function(u){
      vY<-as.numeric(data_mat[,u])
      t(apply(matrix(1:B),1,FUN=function(v){
        p1<-IdxM[[v]]
        Y<-vY[p1]
        Tm<-Tij[p1]
        if(DOP==0){
          Wt<-(1/h)*KernelF(abs(Tm-ETGrid[v])/h)
          if(CovDep){
            cdim<-ifelse(is.matrix(DesignMat),ncol(DesignMat),1)
            if(cdim==1 & !is.matrix(DesignMat)){
              DesignMat<-matrix(DesignMat,ncol=1)
            }
            Xmat<-cbind(rep(1,length(p1)),DesignMat[p1,])
          } else{
            Xmat<-rep(1,length(p1))
          }
          cnp<-as.numeric(sn.mple(x=Xmat,y=Y,w=Wt,penalty = "Qpenalty",opt.method="BFGS")$cp)
          dpp<-as.numeric(cp2dp(cnp,family = "SN"))
          if(CovDep){
            c(cnp[1:(cdim+1)],log(cnp[(cdim+2)]),dpp[(cdim+3)])
          } else{
            c(cnp[1],log(cnp[2]),dpp[3])
          }
        } else{
          if(CovDep){
            cdim<-ifelse(is.matrix(DesignMat),ncol(DesignMat),1)
            if(cdim==1 & !is.matrix(DesignMat)){
              DesignMat<-matrix(DesignMat,ncol=1)
            }
            Xmat<-cbind(rep(1,length(p1)),XMat[p1,v],DesignMat[p1,])
            in_par<-c(mean(Y),runif(1+cdim),log(sd(Y)),runif(1),runif(1))
            as.numeric(optim(par=in_par,fn=LocLogLikSN,t=ETGrid[v],Tij=Tm,Y=Y,meanX=Xmat,scaleX=Xmat,kernel=KernelF,h=h,penalty=PenaltyF)$par)[-c(2,cdim+4)]
          } else{
            Xmat<-cbind(1,XMat[p1,v])
            in_par<-c(mean(Y),runif(1),log(sd(Y)),runif(1),runif(1))
            as.numeric(optim(par=in_par,fn=LocLogLikSN,t=ETGrid[v],Tij=Tm,Y=Y,meanX=Xmat,scaleX=Xmat,kernel=KernelF,h=h,penalty=PenaltyF)$par)[c(1,3,5)]
          }
        }
      }))
    })
    E1datF<-data.frame(sapply(1:length(cds),function(j){
      as.numeric(sapply(locmleres,function(u){u[,j]}))
    }))
    #######
    rm(locmleres)
    ###########
    E1datF$Space<-rep(funARG,each=B)
    E1datF$Time<-rep(ETGrid,times=p)
    fitPLF<-lapply(1:length(cds), function(u){
      frm<-as.formula(paste("X",u,"~",ifelse(cds[u]==2,
                                             "te(Space,Time,k=ES2knots,bs=ES2bs,m=ES2m)",
                                             "s(Space,k=ES2knots[1],bs=ES2bs[1],m=ES2m[1])"),sep=""))
      gam(frm,data = E1datF,method=ES2Esp)
    })
    rm(E1datF)

    ##### Estimation LP process
    LPdat<-data.frame("Space"=rep(funARG,each=length(Tij)),
                      "Time"=rep(Tij,times=p))

    fitPAR<-sapply(seq_len(length(fitPLF)),function(u){
      predict(fitPLF[[u]],newdata=LPdat)
    })

    if(CovDep){
      fullDM<-cbind(1,apply(as.matrix(DesignMat),2,rep,times=p))
    } else{
      fullDM<-1
    }


    lpred<-do.call(cbind,split(rowSums(as.matrix(fitPAR[,1:(length(cds)-2)])*fullDM),LPdat$Space))
    spred<-do.call(cbind,split(exp(fitPAR[,(length(cds)-1)]),LPdat$Space))
    apred<-as.numeric(do.call(cbind,split(fitPAR[,length(cds)],LPdat$Space))[1,])

    rm(LPdat)
    #kapA<-((((1.028571)/h)*2)+1)*(p*B)
    Lh[i]<-sum(sapply(1:p, function(w){
      cpPAR<-t(apply(cbind(lpred[,w],spred[,w],as.numeric(dp2cp(c(0,1,apred[w]),family="SN"))[3]),1,cp2dp,family="SN"))
      dsn(data_mat[,w],xi = cpPAR[,1],omega = cpPAR[,2],alpha = cpPAR[1,3],log = TRUE)
    }))
    if(i ==1){
      i<-i+1
      next
    } else{
      RC<-(abs(Lh[i]-Lh[i-1]))/abs(Lh[i-1])
      i<-i+1
    }
  }
  list("CVresults"=data.frame("Bandwidth"=bandw,"Likelihood"=Lh),
       "OptimalBW"=bandw[length(bandw)],
       "fitOBJ"=fitPLF
       )
}



#' Longitudinal function data generation with skew-normal marginal
#' 
#' Generates longitudinal functional data for n subjects following model introduced in Alam and Staicu (20xx) taking G_alpha to be a Skew-Normal distribution
#'
#' @param argS is a numeric vector contains values of functional argument where the function will be sampled at a given time
#' @param TimePoint is a list of length n gives times points where functional data would have been observed for n subjects
#' @param Sbasis list of length K that has functions generates basis for S as its elements
#' @param Tbasis list of length K where each of the components is a list of length L_k contains basis for time
#' @param Eta is a list of length K where each of the components is a list of length L_k contains eigenvalues
#' @param Sigma2K is a vector of K variance parameters
#' @param Sigma2 is the variance of iid random error
#' @return It returns a list of set of objects
#' \itemize{
#'  \item n: number of subjects generated
#'  \item argS: values of $s_1, s_2,\ldots, s_p$ where the functions were observed
#'  \item Tij: observed time points for the subjects
#'  \item Y: a list of generated data; each element represents data for a subject
#'  \item PWvar: a list of point-wise variance used to generate from U process
#'  \item U: a list of generated U process; elements correspond to subjects
#'  \item W: a list of generated W process; elements correspond to subjects
#' }
#' @example examples/example_data_generation.R
#' @export
SNFData<-function(argS,TimePoint,Sbasis,Tbasis,Eta,Sigma2K,Sigma2,muF,sclF,alpF){
  #require(sn)
  n<-length(TimePoint)
  K<-length(Sbasis)
  Lk<-sapply(Tbasis,length)
  wdata<-lapply(X=1:n, FUN = function(i){
    mi<-length(TimePoint[[i]])
    Reduce(`+`,lapply(1:K, function(k){
      Sb<-Sbasis[[k]](argS)
      Reduce(`+`,lapply(1:(Lk[k]), function(l){
        zeta<-rnorm(1,mean = 0,sd = sqrt(Eta[[k]][l]))
        Tb<-Tbasis[[k]][[l]](TimePoint[[i]])
        zeta*outer(Tb,Sb)
      }))+outer(rnorm(mi,mean = 0,sd = sqrt(Sigma2K[k])),Sb)
    }))+rnorm(mi*length(argS),mean=0,sd=sqrt(Sigma2))
  })
  pwvar<-lapply(X=1:n, FUN = function(i){
    mi<-length(TimePoint[[i]])
    Reduce(`+`,lapply(1:K, function(k){
      Sb<-Sbasis[[k]](argS)
      Reduce(`+`,lapply(1:(Lk[k]), function(l){
        Tb<-Tbasis[[k]][[l]](TimePoint[[i]])
        Eta[[k]][l]*outer(Tb^2,Sb^2)
      }))+Sigma2K[k]*sapply(Sb,function(w){rep(w^2,mi)})
    }))+matrix(rep(Sigma2,mi*length(argS)),nrow=mi,ncol=length(argS))
  })
  udata<-lapply(1:n, function(i){
    mi<-length(TimePoint[[i]])
    tud<-pnorm(wdata[[i]],mean = 0,sd = sqrt(pwvar[[i]]))
    tud[which(tud<0.001,arr.ind = TRUE)]<-0.001
    tud[which(tud>0.999,arr.ind = TRUE)]<-0.999
    tud
  })
  udat<-do.call(rbind,udata)
  subI<-do.call(c,lapply(1:n, function(i){rep(i,length(TimePoint[[i]]))}))
  allTP<-do.call(c,TimePoint)
  ydata<-split.data.frame(do.call(cbind,lapply(X=1:length(argS), FUN=function(l){
    s<-argS[l]
    msPar<-t(sapply(allTP,function(t){c(muF(s,t),sclF(s,t))}))
    msPar[,1]+(msPar[,2]*qsn(udat[,l],dp=cp2dp(c(0,1,as.numeric(dp2cp(c(0,1,alpF(s)),family = "SN"))[3]),family = "SN"),solver = "RFB"))
  })),subI)
  list("n"= n,
       "argS"=argS,
       "Tij"=TimePoint,
       "Y"=ydata,
       "PWvar"=pwvar,
       "U"=udata,
       "W"=wdata)
}


#' Longitudinal functional data analysis using FPCA
#' 
#' Implements longitudinal functional data analysis (Park and Staicu, 2015). It decomposes longitudinally-observed functional observations in two steps. It first applies FPCA on a properly defined marginal covariance function and obtain estimated scores (mFPCA step). Then it further models the underlying process dynamics by applying another FPCA on a covariance of the estimated scores obtained in the mFPCA step. The function also allows to use a random effects model to study the underlying process dynamics instead of a KL expansion model in the second step. Scores in mFPCA step are estimated using numerical integration. Scores in sFPCA step are estimated under a mixed model framework.
#'
#' @param Y	a matrix of which each row corresponds to one curve observed on a regular and dense grid (dimension of N by m; N = total number of observed functions; m = number of grid points)
#' @param subject.index	subject id; vector of length N with each element corresponding a row of Y
#' @param visit.index	index for visits (repeated measures); vector of length N with each element corresponding a row of Y
#' @param obsT	actual time of visits at which a function is observed; vector of length N with each element corresponding a row of Y
#' @param funcArg	numeric; function argument
#' @param numTEvalPoints	total number of equidistant time points where mean function is estimated, defaults to 41
#' @param newdata	an optional data frame providing predictors (i for subject id / Ltime for visit time) with which prediction is desired; defaults to NULL
#' @param fbps.knots	list of two vectors of knots or number of equidistanct knots for all dimensions for a fast bivariate P-spline smoothing (fbps) method used to estimate a bivariate, smooth mean function; defaults to c(5,10); see fbps
#' @param fbps.p	integer;degrees of B-spline functions to use for a fbps method; defaults to 3; see fbps
#' @param fbps.m	integer;order of differencing penalty to use for a fbps method; defaults to 2; see fbps
#' @param mFPCA.pve	proportion of variance explained for a mFPCA step; used to choose the number of principal components (PCs); defaults to 0.95; see fpca.face
#' @param mFPCA.knots	number of knots to use or the vectors of knots in a mFPCA step; used for obtain a smooth estimate of a covariance function; defaults to 35; see fpca.face
#' @param mFPCA.p	integer; the degree of B-spline functions to use in a mFPCA step; defaults to 3; see fpca.face
#' @param mFPCA.m	integer;order of differencing penalty to use in a mFPCA step; defaults to 2; see fpca.face
#' @param mFPCA.npc	pre-specified value for the number of principal components; if given, it overrides pve; defaults to NULL; see fpca.face
#' @param LongiModel.method	model and estimation method for estimating covariance of estimated scores from a mFPCA step; either KL expansion model or random effects model; defaults to fpca.sc
#' @param sFPCA.pve proportion of variance explained for sFPCA step; used to choose the number of principal components; defaults to 0.95; see fpca.sc
#' @param sFPCA.nbasis	number of B-spline basis functions used in sFPCA step for estimation of the mean function and bivariate smoothing of the covariance surface; defaults to 10; see fpca.sc
#' @param sFPCA.npc	pre-specified value for the number of principal components; if given, it overrides pve; defaults to NULL; see fpca.sc
#' @param gam.method	smoothing parameter estimation method when gam is used for predicting score functions at unobserved visit time, T; defaults to REML; see gam
#' @param gam.kT	dimension of basis functions to use; see gam
#' @return a list with the following items:
#' \itemize{
#'  \item obsData: a list of	observed data (input), i: subject id, j: visit index, Tij: observed time points, funcArg: functional argument
#'  \item funcArg: function argument
#'  \item visitTime	: visit times
#'  \item fitted.values	: fitted values (in-sample); of the same dimension as Y
#'  \item fitted.values.all	: a list of which each component consists of a subject's fitted values at all pairs of evaluation points (s and T)
#'  \item predicted.values : predicted values for variables provided in newdata
#'  \item bivariateSmoothMeanFunc	: estimated bivariate smooth mean function
#'  \item meanFOBJ : fitted object for bi-variate mean function
#'  \item mcMEAN : fitted mean for marginal covariance
#'  \item mFPCA.efunctions : estimated eigenfunction in a mFPCA step
#'  \item mFPCA.evalues	: estimated eigenvalues in a mFPCA step
#'  \item mFPCA.npc	: number of principal components selected with pre-specified pve in a mFPCA step
#'  \item mFPCA.scree.eval : estimated eigenvalues obtained with pre-specified pve = 0.9999; for scree plot
#'  \item sFPCA.xiHat.bySubj :	a list of which each component consists of a subject's predicted score functions evaluated at equidistanced grid in direction of visit time, T
#'  \item sFPCA.npc	: a vector of numbers of principal components selected in a sFPCA step with pre-specified pve; length of mFPCA.npc
#'  \item mFPCA.covar	: estimated marginal covariance
#'  \item sFPCA.longDynCov.k : a list of estimated covariance of score function; length of mFPCA.npc
#'  \item sigma2k : estimated variance for every k
#' }
#' @example examples/example_lfda.R
#' @export
fpcaLFDA<-function (Y, subject.index, visit.index, obsT = NULL, funcArg = NULL,
                    numTEvalPoints = 41, newdata = NULL, fbps.knots = c(5, 10),
                    fbps.p = 3, fbps.m = 2, mFPCA.pve = 0.95, mFPCA.knots = 35,
                    mFPCA.p = 3, mFPCA.m = 2, mFPCA.npc = NULL, LongiModel.method = c("fpca.sc", "lme"), sFPCA.pve = 0.95, sFPCA.nbasis = 10,
                    gam.method = "REML", gam.kT = 10)
{
  y <- as.matrix(Y)
  colnames(y) <- NULL
  if (is.null(funcArg)) {
    ss <- seq(0, 1, length.out = ncol(y))
  }else {
    if (length(funcArg) == ncol(y)) {
      ss <- funcArg
    }else {
      warning("length(funcArg) and ncol(Y) do not match. funcArg is re-defined: funcArg = ncol(Y) equally spaced grid points in [0,1].")
      ss <- seq(0, 1, length.out = ncol(y))
    }
  }
  Tij <- obsT
  ETGrid <- seq(min(Tij), max(Tij), length.out = numTEvalPoints)
  TT <- sort(unique(Tij))
  numTEvalPoints<-length(TT)
  n <- length(unique(subject.index))
  mi <- aggregate(subject.index, by = list(subject.index),
                  length)[, 2]
  subject.index <- unlist(sapply(1:n, function(a) rep(a, mi[a])))
  M <- length(ss)
  J <- length(TT)
  Ncurves <- nrow(y)
  uTij <- unique(Tij)
  fit.fbps <- fbps(data = y, covariates = list(Tij, ss), knots = fbps.knots)
  mu.hat <- fit.fbps$Yhat
  new.y <- y - mu.hat
  fit1 <- fpca.face(Y = new.y,argvals = ss, pve = mFPCA.pve, knots = mFPCA.knots,
                    p = mFPCA.p, m = mFPCA.m, npc = mFPCA.npc)
  phi.hat <- fit1$efunctions * sqrt(M)
  K.hat <- fit1$npc
  lambda.hat <- fit1$evalues/M
  marCovar.hat <- Reduce("+", lapply(seq_len(length(lambda.hat)),
                                     function(a) lambda.hat[a] * t(t(phi.hat[, a])) %*% t(phi.hat[,
                                                                                                  a])))
  fitfit <- fpca.face(Y = new.y, pve = 0.9999)
  scree.PVE <- cumsum(fitfit$evalues)/sum(fitfit$evalues)
  v <- unlist(lapply(Tij, function(a) which(abs(TT - a) ==
                                              min(abs(TT - a)))))
  if (LongiModel.method == "fpca.sc") {
    xi.hat0 <- list()
    for (k in seq_len(K.hat)) {
      xi.hat0[[k]] <-data.frame(".id"=subject.index,".index"=obsT,".value"=fit1$scores[, k]/sqrt(M))
    }
    fit2 <- lapply(xi.hat0, function(a) fpca.sc(ydata = a, pve = sFPCA.pve,
                                                var = TRUE))
    sigma2k<-sapply(fit2,function(u){u$sigma2})
    xi.hat <- lapply(fit2, function(a) a$Yhat)
    longDynamicsCov.hat.k <- lapply(seq_len(fit1$npc), function(a) Reduce("+",
                                                                          lapply(seq_len(length(fit2[[a]]$evalues)), function(b) fit2[[a]]$evalues[b] *
                                                                                   t(t(fit2[[a]]$efunctions[, b])) %*% t(fit2[[a]]$efunctions[,
                                                                                                                                              b]))))
  }
  else if (LongiModel.method == "lme") {
    lme.coef <- lme.pred <- lme.full.pred <- lme.fit <- list()
    sigma2k<-NULL
    for (k in seq_len(K.hat)) {
      lme.dat <- data.frame(Y = fit1$scores[, k]/sqrt(M),
                            X = Tij, subj = subject.index)
      lme.fit0 <- lmer(Y ~ (X | subj), data = lme.dat,
                       REML = TRUE)
      sigma2k[k]<-summary(lme.fit0)$sigma^2
      lme.fit[[k]] <- lme.fit0
      lme.coef[[k]] <- coef(lme.fit0)
      lme.pred[[k]] <- fitted(lme.fit0)
      lme.full.pred[[k]] <- t(matrix(predict(lme.fit0,
                                             newdata = data.frame(X = rep(TT, n), subj = rep(unique(subject.index),
                                                                                             each = length(TT)))), length(TT)))
    }
    xi.hat <- lme.full.pred
    longDynamicsCov.hat.k <- lapply(seq_len(fit1$npc), function(a) cbind(rep(1,
                                                                             numTEvalPoints), TT) %*% as.matrix(as.data.frame(VarCorr(lme.fit[[a]])[[1]])) %*%
                                      t(cbind(rep(1, numTEvalPoints), TT)))
  }
  fitted <- mu.hat + t(matrix(rep(fit1$mu, Ncurves), nrow = length(ss))) +
    do.call(rbind, lapply(seq_len(Ncurves), function(icv) unlist(lapply(xi.hat,
                                                                        function(a) a[subject.index[icv], v[icv]])))) %*%
    t(phi.hat)
  muHat <- matrix(predict(object = fit.fbps, newdata = data.frame(x = rep(ETGrid,
                                                                          M), z = rep(ss, each = length(ETGrid))))$fitted.values,
                  nrow = length(ETGrid))

  xi.hat.bySubj <- lapply(1:n, function(i) sapply(xi.hat, function(a) a[i,
  ]))
  xi.hat.phi.hat.bySubj <- lapply(xi.hat.bySubj, function(a){
    Amat<-apply(a, 2,FUN = function(v){approx(x=TT,y = v,xout = ETGrid,rule = 2)$y})
    t(apply(Amat,1, function(b) matrix(b, nrow = 1) %*% t(phi.hat)))})
  Yhat.all <- lapply(xi.hat.phi.hat.bySubj, function(a) a +
                       muHat + t(matrix(rep(fit1$mu, length(ETGrid)), nrow = length(ss))))
  if (!is.null(newdata)) {
    Jpred <- nrow(newdata)
    randDev <- function(row) {
      i <- newdata$i[row]
      Tpred <- newdata$Ltime[row]
      temp.xi.hat.bySubj <- xi.hat.bySubj[[i]]
      temp.fit <- apply(xi.hat.bySubj[[i]], 2, function(a) {
        temp.data <- data.frame(y = a, x = TT)
        gam(y ~ s(x, k = gam.kT, bs = "cr"), data = temp.data,
            method = gam.method)
      })
      xi.hat.atT <- sapply(temp.fit, function(a) predict.gam(a,
                                                             newdata = data.frame(x = Tpred)))
      return(xi.hat.atT %*% t(phi.hat))
    }
    muHat.pred <- matrix(predict.fbps(object = fit.fbps,
                                      newdata = data.frame(x = rep(newdata$Ltime, M), z = rep(ss,
                                                                                              each = Jpred)))$fitted.values, nrow = Jpred)
    predicted <- muHat.pred + t(matrix(rep(fit1$mu, Jpred),
                                       nrow = length(ss))) + do.call(rbind, lapply(seq_len(Jpred),
                                                                                   function(icv) randDev(icv)))
  }
  else {
    predicted <- NULL
    newdata.hat <- NULL
  }
  if (LongiModel.method == "fpca.sc") {
    ret <- list(obsData = list(y = Y, i = subject.index,
                               j = visit.index, Tij = obsT, funcArg = ss), i = subject.index,
                funcArg = ss, visitTime = TT, fitted.values = fitted,
                fitted.values.all = Yhat.all, predicted.values = predicted,
                bivariateSmoothMeanFunc = muHat,meanFOBJ=fit.fbps, mcMEAN=fit1$mu, mFPCA.efunctions = phi.hat,
                mFPCA.evalues = lambda.hat, mFPCA.npc = K.hat, mFPCA.scree.eval = fitfit$evalues/M,
                sFPCA.xiHat.bySubj = xi.hat.bySubj, sFPCA.npc = unlist(lapply(fit2,
                                                                              function(a) a$npc)), mFPCA.covar = marCovar.hat,
                sFPCA.longDynCov.k = longDynamicsCov.hat.k,
                sigma2k=sigma2k)
  }
  else if (LongiModel.method == "lme") {
    ret <- list(obsData = list(y = Y, i = subject.index,
                               j = visit.index, Tij = obsT, funcArg = ss), i = subject.index,
                funcArg = ss, visitTime = TT, fitted.values = fitted,
                fitted.values.all = Yhat.all, predicted.values = predicted,
                bivariateSmoothMeanFunc = muHat, meanFOBJ=fit.fbps, mcMEAN=fit1$mu, mFPCA.efunctions = phi.hat,
                mFPCA.evalues = lambda.hat, mFPCA.npc = K.hat, mFPCA.scree.eval = fitfit$evalues/M,
                sFPCA.xiHat.bySubj = xi.hat.bySubj, mFPCA.covar = marCovar.hat,
                sFPCA.longDynCov.k = longDynamicsCov.hat.k,
                sigma2k=sigma2k)
  }
  class(ret) <- "lfpca"
  return(ret)
}

#' In-sample prediction via LFDA
#' 
#' Perform prediction using the results from fpcaLFDA function; prediction is for subjects in the training sample
#'
#' @param lfdaOBJ fitted object fpcaLFDA
#' @param TimeOBJ is a vector of observed time points stacked for all subjects
#' @param gridS is a vector of functional argument where trajectories will be predicted
#' @param gTID is the subject id for which the prediction will be made
#' @param gridT is a list of length equals length of gTID contains time points where prediction will be made for the subjects identified in gTID
#' @return a list with following items
#' \itemize{
#'  \item PredFD : a list with predicted function for subjects
#'  \item Subjects : id of subjects (gTID)
#'  \item funARG : argument where trajectories were predicted (gridS)
#'  \item TimePoints : a list of time where prediction is (gridT)
#'  \item Sbasis : basis function estimated from marginal covariance
#'  \item TbasisOBJ : basis function for time-varying coefficients
#' }
#' @example examples/example_lfda_ISpred.R
#' @export
PWPRED<-function(lfdaOBJ,TimeOBJ,gridS=NULL,gTID=NULL,gridT=NULL){
  xiHAT<-lfdaOBJ$sFPCA.xiHat.bySubj
  Sbasis<-lfdaOBJ$mFPCA.efunctions
  mOBJ<-lfdaOBJ$meanFOBJ
  argS<-lfdaOBJ$funcArg
  sc1mu<-approxfun(x=argS,y=lfdaOBJ$mcMEAN,rule = 2)
  if(!is.null(gridS))
    Sbasis<-apply(Sbasis, 2, FUN=function(u){
      approx(x=argS,y=u,xout = gridS,rule = 2)$y
    })
  if(is.null(gridS))
    gridS<-argS
  n<-length(xiHAT)
  K<-ncol(Sbasis)
  predSCR<-do.call(what=rbind,args=xiHAT)
  XiHat<-lapply(as.list(1:K),FUN=function(u){matrix(predSCR[,u],ncol=length(unique(TimeOBJ)),byrow=TRUE)})
  PsiH<-lapply(as.list(1:K),FUN=function(u){
    fpca.sc(Y=XiHat[[u]],argvals = sort(unique(TimeOBJ)))
  })

  # Mean function for time-varying coefficients
  psiMU<-lapply(PsiH,function(u){
    approxfun(x=u$argvals,u$mu,rule = 2)
  })

  if(is.null(gTID)){
    gTID<-1:n
    PredFD<-lapply(seq_len(length(gTID)), function(i){
      Reduce(`+`,lapply(seq_len(K), function(k){
        Lk<-ncol(PsiH[[k]]$efunctions)
        Reduce(`+`,lapply(seq_len(Lk),FUN=function(l){
          (PsiH[[k]]$scores[gTID[i],l])*outer(approx(x=sort(unique(TimeOBJ)),y=PsiH[[k]]$efunctions[,l],xout = gridT[[i]],rule = 2)$y,Sbasis[,k])
        })) + outer(psiMU[[k]](gridT[[i]]),Sbasis[,k])
      })) + outer(rep(1,length(gridT[[i]])),sc1mu(gridS)) + matrix(predict(object = mOBJ,
                                newdata = data.frame(x = rep(gridT[[i]],length(argS)),
                                                     z = rep(argS, each = length(gridT[[i]]))))$fitted.values,
                   nrow = length(gridT[[i]]))
    })
  } else{
    PredFD<-lapply(seq_len(length(gTID)), function(i){
      Reduce(`+`,lapply(seq_len(K), function(k){
        Lk<-ncol(PsiH[[k]]$efunctions)
        Reduce(`+`,lapply(seq_len(Lk),FUN=function(l){
          (PsiH[[k]]$scores[gTID[i],l])*outer(approx(x=sort(unique(TimeOBJ)),y=PsiH[[k]]$efunctions[,l],xout = gridT[[i]],rule = 2)$y,Sbasis[,k])
        })) + outer(psiMU[[k]](gridT[[i]]),Sbasis[,k])
      })) + outer(rep(1,length(gridT[[i]])),sc1mu(gridS)) +matrix(predict(object = mOBJ,
                              newdata = data.frame(x = rep(gridT[[i]],length(argS)),
                                                   z = rep(argS, each = length(gridT[[i]]))))$fitted.values,
                 nrow = length(gridT[[i]]))
    })
  }
  # Return object
  list("PredFD"=PredFD,
       "Subjects"=gTID,
       "funARG" = gridS,
       "TimePoints"=gridT,
       "Sbasis"=Sbasis,
       "TbasisOBJ"=PsiH)
}

#' Out-of-sample prediction via the LFDA
#' 
#' Perform prediction using the results from fpcaLFDA function
#'
#' @param lfdaOBJ fitted object fpcaLFDA
#' @param TimeOBJ is a vector of observed time points stacked for all subjects
#' @param nfunDATA a list of functional data observed for new subjects
#' @param nfunARG argument where functions for new subjects were observed
#' @param nobsTIME a list of time points; every element represents the observed time points for subjects in the nfunDATA
#' @param identical.ARG a logical scalar indicates whether observed functional grid for new subjects is same as the subjects used in the training the LFDA
#' @param gridS is a vector of functional argument where trajectories will be predicted
#' @param gridT is a list of length equals length of gTID contains time points where prediction will be made for the subjects identified in gTID
#' @return a list with following items
#' \itemize{
#'  \item PredFD : a list with predicted function for subjects
#'  \item Subjects : id of subjects (gTID)
#'  \item funARG : argument where trajectories were predicted (gridS)
#'  \item TimePoints : a list of time where prediction is (gridT)
#'  \item Sbasis : basis function estimated from marginal covariance
#'  \item TbasisOBJ : basis function for time-varying coefficients
#' }
#' @example examples/example_lfda_OSpred.R
#' @export
OutSamPWPRED<-function(lfdaOBJ,TimeOBJ,nfunDATA,nfunARG=NULL,nobsTIME,
                       identical.ARG=TRUE,gridS=NULL,gridT){
  xiHAT<-lfdaOBJ$sFPCA.xiHat.bySubj
  phiH<-lfdaOBJ$mFPCA.efunctions
  mOBJ<-lfdaOBJ$meanFOBJ
  argS<-lfdaOBJ$funcArg
  nmi<-sapply(nobsTIME,length)
  nID<-rep(1:length(nmi),nmi)
  sc1mu<-approxfun(x=argS,y=lfdaOBJ$mcMEAN,rule = 2)
  if(!identical.ARG){
    phiH<-apply(phiH, 2, FUN = function(u){
      approx(x=funARG,y=u,xout = nfunARG,rule=2)$y
    })
    if(is.null(nfunARG))
      stop(print("For non-identical functional argument, provide the value of nfunARG"))
    argS<-nfunARG
  } else{
    if(is.null(nfunARG))
      nfunARG<-argS
  }

  if(!is.null(gridS))
    phiH<-apply(phiH, 2, FUN=function(u){
      approx(x=argS,y=u,xout = gridS,rule = 2)$y
    })
  n<-length(xiHAT)
  nP<-length(nfunDATA)
  Khat<-ncol(phiH)
  predSCR<-do.call(what=rbind,args=xiHAT)
  XiHat<-lapply(as.list(1:Khat),FUN=function(u){matrix(predSCR[,u],ncol=length(unique(TimeOBJ)),byrow=TRUE)})
  PsiH<-lapply(as.list(1:Khat),FUN=function(u){
    fpca.sc(Y=XiHat[[u]],argvals = sort(unique(TimeOBJ)),var = TRUE)
  })
  psiKL<-lapply(1:Khat, function(k){
    Lk<-PsiH[[k]]$npc
    lapply(1:Lk, function(l){
      approxfun(x=PsiH[[k]]$argvals,y=PsiH[[k]]$efunctions[,l],rule = 2)
    })
  })

  psiMU<-lapply(PsiH,function(u){
    approxfun(x=u$argvals,u$mu,rule = 2)
  })

  #############Score estimation for new subjects################
  dSW<-lapply(seq_len(length(nfunDATA)), function(u){
    nfunDATA[[u]]-matrix(predict(object = mOBJ,
                                      newdata = data.frame(x = rep(nobsTIME[[u]],length(argS)),
                                                           z = rep(argS, each = nmi[u])))$fitted.values,
                         nrow = nmi[u])
  })


  NXi<-do.call(rbind,lapply(seq_len(length(dSW)), function(w){
    if(length(nobsTIME[[w]])==1){
      apply(phiH,2,FUN = function(v){
        sum((dSW[[w]]-sc1mu(nfunARG))*v*diff(c(0,nfunARG)))
      })
    } else{
      t(sapply(seq_len(length(nobsTIME[[w]])), function(j){
        apply(phiH,2,FUN = function(v){
          sum((dSW[[w]][j,]-sc1mu(nfunARG))*v*diff(c(0,nfunARG)))
        })
      }))
    }
  }))
  subGK<-lapply(seq_len(length(nobsTIME)),function(i){
    lapply(1:Khat,function(k){
      Lk<-PsiH[[k]]$npc
      Gk<-Reduce(`+`,lapply(1:Lk, function(l){
        PsiH[[k]]$evalues[l]*outer(psiKL[[k]][[l]](nobsTIME[[i]]),psiKL[[k]][[l]](nobsTIME[[i]]))
      }))
      if(lfdaOBJ$sigma2k[k]<1e-10){
        solve(Gk+(1e-10*diag(length(nobsTIME[[i]]))))
      } else{
        solve(Gk+(lfdaOBJ$sigma2k[k]*diag(length(nobsTIME[[i]]))))
      }
    })
  })
  nZETA<-lapply(seq_len(Khat),function(k){
    knXI<-split(NXi[,k],nID)
    Lk<-PsiH[[k]]$npc
    do.call(rbind,lapply(seq_len(length(nobsTIME)),function(i){
      sapply(1:Lk, function(l){
        PsiH[[k]]$evalues[l]*sum(psiKL[[k]][[l]](nobsTIME[[i]])*(subGK[[i]][[k]]%*%matrix(knXI[[i]]-psiMU[[k]](nobsTIME[[i]]),ncol=1)))
      })
    }))
  })

  PredFD<-lapply(seq_len(nP), function(i){
    Reduce(`+`,lapply(seq_len(Khat), function(k){
      Lk<-PsiH[[k]]$npc
      Reduce(`+`,lapply(seq_len(Lk),FUN=function(l){
        (nZETA[[k]][i,l])*outer(approx(x=sort(unique(TimeOBJ)),y=PsiH[[k]]$efunctions[,l],xout = gridT[[i]],rule = 2)$y,phiH[,k])
      })) + outer(psiMU[[k]](gridT[[i]]),phiH[,k])
    })) + outer(rep(1,length(gridT[[i]])),sc1mu(nfunARG)) + matrix(predict(object = mOBJ,
                              newdata = data.frame(x = rep(gridT[[i]],length(argS)),
                                                   z = rep(argS, each = length(gridT[[i]]))))$fitted.values,
                 nrow = length(gridT[[i]]))
  })
  list("PredFD"=PredFD,
       "Subjects"=unique(nID),
       "funARG" = nfunARG,
       "TimePoints"=gridT,
       "Sbasis"=phiH,
       "TbasisOBJ"=PsiH)
}


#' Epanechnikov kernel
#' 
#' Epanechnikov Kernel function with band-width 1;
#' for a bandwidth of h, the input should be scaled by h
#'
#' @param x point at which kernel will be evaluated
#' @returns a scalar of kernel value
#' @examples
#' depan(0.3)
#'
#' @export
depan<-function(x){
  ifelse(abs(x)<=1,(3/4)*(1-x^2),0)
}

