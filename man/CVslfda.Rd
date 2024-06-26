% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SkewedFDA.R
\name{CVslfda}
\alias{CVslfda}
\title{Cross-validation in two-step estimation of PLF by the sLFDA}
\usage{
CVslfda(
  funDATA,
  funARG,
  obsTIME,
  ETGrid,
  DOP = 1,
  KernelF,
  Hgrid,
  PenaltyF = Qpenalty,
  ES2knots = c(10, 10),
  ES2bs = c("ps", "ps"),
  ES2m = c(2, 2),
  ES2Esp = "REML",
  Thresh = 0.05,
  CovDep = FALSE,
  DesignMat = NULL,
  CovDepStr = NULL,
  parCOMP = FALSE,
  n_cores = 7,
  par_seed = 100
)
}
\arguments{
\item{funDATA}{is a list of functional data observed from different subjects}

\item{funARG}{functional arguments where the underlying function}

\item{obsTIME}{is a list of time points for subjects in the study}

\item{ETGrid}{is a vector of regularly spaced time points}

\item{DOP}{is a scalar represents degree of polynomial for local fitting for mean and scale function}

\item{KernelF}{is a character represents the name of the standard kernel to be used local linear maximum likelihood estimation}

\item{Hgrid}{a numeric vector of bandwidth values for which log-likelihood will be evaluated}

\item{PenaltyF}{is the penalty function for shape parameter of the skewnormal distribution}

\item{ES2knots}{is a  numeric vector of length two represent the numbers of knots for functional and longitudinal domain required for smoothing at step 2 of the estimation}

\item{ES2bs}{is a vector of two strings represent names of basis functions in S and T}

\item{ES2m}{order of penalities to be used in step of PLF estimation}

\item{ES2Esp}{is a string represents the name of method to be used in step 2 for determining tuning parameter}

\item{Thresh}{is a scalar represents the minimum increment for likelihood with change in bandwidth}

\item{CovDep}{is a logical scalar that represents whether mean depends on covariates or not}

\item{DesignMat}{is a matrix of baseline covariate for subjects in the data}

\item{CovDepStr}{is a vecotr of length equals to the number of covariates.
Each element can be either 1 or 2; 1 represents coefficient depends on s only, 2 means
the coefficient depends on both s and t}

\item{parCOMP}{a logical scalar indicates whether computation is performed parallely or not}

\item{n_cores}{number of cores to use while performing the computation in parallel}

\item{par_seed}{a scalar used as seed in ``dorng'' when parallel computing would be used; for sequential computation, skewedFDA should preceed set.seed for reproducibility}

\item{CVThresh}{is a scalar represents the minimum increment for likelihood with change in bandwidth}
}
\value{
a list with the following items:
\itemize{
 \item CVresults : a data frame with two columns (bandwidth and log-likelihood)
 \item OptimalBW : chosen optimal bandwidth value
 \item fitOBJ : smooth fitted objects for population level functions at the optimal bandwidth
 \item AvgTime : Mean time (in seconds) taken to obtain PLF for a given kernel bandwidth
 \item TotalTime : Total time (in seconds) taken to obtain PLF in the cross-validation
}
}
\description{
Cross-Validation based on likelihood function for bandwidth selection at the
first step of estimation proposed for population level function estimation in
Alam and Staicu (20xx)
}
\examples{
# This an example of skewed longitudinal functional data where marginal
# is a skew-normal distribution
# Required function for Fourier basis and skew-normal distribution
library(fda)
library(sn)
library(mgcv)
# Population level functions
# Mean function
meanPF<-function(s,t){
  2.5+(3.15*s)+(4*t)+(2*s*t)
}

# Scale Function
sFUN<-function(s,t){
  require(mvtnorm)
  25*dmvnorm(c(s,t),mean = c(0,0),sigma=matrix(c(2.5,0.75,0.75,3.5),2,2))
}

# Shape Function

alFUN<-function(s){
  2e1*(exp(2*s)/(1+exp(2*s)))*sin(6*pi*s/4)
}

# Sbasis functions for $L^2[S]$
SFbasis<-list("sfourier1"=function(s){
  fourier(x=s,nbasis=5,period = 1)[,2]
},"sfourier2"=function(s){
  fourier(x=s,nbasis=5,period = 1)[,3]
},"sfourier3"=function(s){
  fourier(x=s,nbasis=5,period = 1)[,4]
})

# Basis functions for time-dependent coefficients

TBasis<-list(list("tb11" = function(T){fourier(x=T,nbasis=3,period=1)[,2]},
                  "tb12" = function(T){fourier(x=T,nbasis=3,period=1)[,3]}),
             list("tb21" = function(T){fourier(x=T,nbasis=5,period=1)[,4]},
                  "tb22" = function(T){fourier(x=T,nbasis=5,period=1)[,5]}),
             list("tb31" = function(T){fourier(x=T,nbasis=7,period=1)[,6]},
                  "tb32" = function(T){fourier(x=T,nbasis=7,period=1)[,7]})
)


# A set of regular grid in functional domain S
n_s<-51
ss<-seq(0,1,length.out = n_s)

# Eigen values for time-dependent coefficient $eta_{kl}$
tdcfV<-list(c(0.6,0.4),c(0.5,0.3),c(0.25,0.20))


# Variance parameters
sig2_S<-c(0.328,0.210,0.046)
sig2_E<-0.183


# Testing 
ntp<-20
tp<-matrix(seq(0,1,length.out = ntp+2)[-c(1,ntp+2)],ncol=1)

# Full process time grid
Tg<-seq(0,1,length.out=51)

## number of subjects
n<-25

# Selection of mi
mi<-sapply(1:n,sample,x=6:12,size = 1)

# Observed Time points
Tij<-lapply(seq_len(length(mi)), function(w){sort(sample(Tg,mi[w]))})

# Simulation of Skewed FD
gdata<-SNFData(argS = ss,TimePoint = Tij,Sbasis = SFbasis,Tbasis = TBasis,
               Eta = tdcfV,Sigma2K = sig2_S,Sigma2 = sig2_E,
               muF = meanPF,sclF = sFUN,alpF = alFUN)


# cross-validation to choose a bandwidth at step 1

cv_fit<-CVslfda(funDATA=gdata$Y,funARG=ss,obsTIME=Tij,ETGrid=tp,DOP=1,KernelF=depan,
                Hgrid=seq(0.02,0.2,0.02),PenaltyF=Qpenalty,ES2knots=c(10,10),
                ES2bs=c("ps","ps"),ES2m=c(2,2),ES2Esp="REML",Thresh=0.05,
                CovDep=FALSE,DesignMat=NULL,CovDepStr=NULL)

}
