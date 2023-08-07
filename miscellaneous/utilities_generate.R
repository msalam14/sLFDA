# Generates longitudinal functional data for n subjects following model introduced in Alam and Staicu (20xx) taking G_alpha to be a Skew-Normal distribution
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
#' @example example_data_generation.R
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