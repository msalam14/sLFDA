# This an example of skewed longitudinal functional data where marginal
# is a skew-normal distribution
# Required function for Fourier basis and skew-normal distribution
library(fda)
library(sn)
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
n<-2

# Selection of mi
mi<-sapply(1:n,sample,x=3:6,size = 1)

# Observed Time points
Tij<-lapply(seq_len(length(mi)), function(w){sort(sample(Tg,mi[w]))})

# Simulation of Skewed FD
gdata<-SNFData(argS = ss,TimePoint = Tij,Sbasis = SFbasis,Tbasis = TBasis,
               Eta = tdcfV,Sigma2K = sig2_S,Sigma2 = sig2_E,
               muF = meanPF,sclF = sFUN,alpF = alFUN)
