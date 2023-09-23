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

# Bi-variate basis
STBasis<-list("stB1"=function(s,t){
  sapply(s, function(u){
    0.8578*sin(u^2+(0.5*t^2))
  })
},
"stB2"=function(s,t){
  sapply(s, function(u){
    (0.8721*sin((0.3*u^2)+(0.6*t^2)))-(0.2988*sin(u^2+(0.5*t^2)))
  })
})

# Marginal basis

SBasis<-list("sB1"=function(s){
  0.8578*sin(s^2)
},
"sB2"=function(s,t){
  (0.8721*sin(0.3*s^2))-(0.2988*sin(s^2))
})




# A set of regular grid in functional domain S
n_s<-51
ss<-seq(0,1,length.out = n_s)

# Eigen values for time-dependent coefficient $eta_{kl}$
tdcfV<-c(1.6,0.8)


# Variance parameters
sig2_S<-c(0.025,0.015)
sig2_E<-0.08


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
gdata<-SNFDataBB(argS = ss,TimePoint = Tij,STbasis = STBasis,Sbasis = SBasis,
               Eta = tdcfV,Sigma2K = sig2_S,Sigma2 = sig2_E,
               muF = meanPF,sclF = sFUN,alpF = alFUN)
