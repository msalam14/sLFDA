# This is an example of log-likelihood using local linear approximation
library(sn)
set.seed(125)
Y<-rsn(n=50,xi = 5.5,omega = 1,alpha = 0.3)
Tij<-round(sort(runif(50)),2)
t<-0.25
Xmat<-cbind(rep(1,50),(t-Tij))


LocLogLikSN(runif(5),t=t,Tij=Tij,Y=Y,meanX=Xmat,scaleX=Xmat,kernel=depan,h=0.12,penalty=Qpenalty)