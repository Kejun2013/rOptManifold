#####Two-step Methods#####
#####generate basis####
library(splines)
#two principal components
set.seed(100)
vv=function(t){ 
  f1=sin(pi*t)/sqrt(5)
  f2=cos(pi*t)/sqrt(5)
  cbind(f1,f2)
}
x=seq(0,1,by=0.001)
truFun=vv(x)*sqrt(10)

#t(truFun)%*%truFun*0.001
#plot(truFun[,2],type="l")
#Basis
M=10
l=1/(M+1)
knots=c(c(-2*l,-l),seq(0,1,length.out=M+2),c(1+l,1+2*l))
D=splineDesign(knots=knots,x=x,ord=3)
D2=splineDesign(knots=knots,x=x,ord=3,derivs=rep(2,length(x)))
#####orthorgonalization####
Basis=qr.Q(qr(D))*sqrt(1000)
dim(Basis)
R=qr.R(qr(D))
BasisDeriv=D2%*%solve(R)
Omega=t(BasisDeriv)%*%BasisDeriv/1000^2
n=dim(Omega)[1]

##generate samples
nSample=300 #Sample Size
ySample=list();tSample=list();bSample=list();
for(i in 1:nSample){
  ni=sample(10:20,1)
  tTemp=sample(0:1000,ni,replace=FALSE)
  bSample=c(bSample,list(Basis[tTemp+1,]/ni))
  tTemp=tTemp/1000
  tSample=c(tSample,list(tTemp))
  yTemp=vv(tTemp)%*%matrix(c(6*rnorm(1),3*rnorm(1)),ncol=1)+
              matrix(rnorm(ni,sd=1),ncol=1)
  ySample=c(ySample,list(yTemp/ni) )
}


objF=function(Sigma){
  resiSum=0
  for(ii in 1:nSample){
    temp=ySample[[ii]]%*%t(ySample[[ii]])-bSample[[ii]]%*%Sigma%*%t(bSample[[ii]])
    #if(ii==10) print(temp[1:3,1:3])
    resiSum=resiSum+sum(temp^2)
  }
  objValue=resiSum+lambda*sum(Omega*Sigma)
  objValue
}

gradF=function(Sigma){
#     Sigma=matrix(runif(2*n),n,2)
#     Sigma=Sigma%*%t(Sigma)
  
  temp=matrix(0,n,n)
  for(ii in 1:nSample){
    temp=temp+t(bSample[[ii]])%*%(-ySample[[ii]]%*%t(ySample[[ii]])+
                      bSample[[ii]]%*%Sigma%*%t(bSample[[ii]]))%*%bSample[[ii]]
  }
  temp=-2*temp+Omega*lambda
  #objF(Sigma-0.001*temp)
  temp
}

hessF=function(Sigma,Z){
  #Sigma=problem@Y[[1]]#
  temp=matrix(0,n,n)
  for(ii in 1:nSample){
    BtB=t(bSample[[ii]])%*%bSample[[ii]]
    temp=temp+BtB%*%Z%*%BtB
  }
  temp=2*temp
  #objF(Sigma-0.001*temp)
  temp
}





problem=fixedRankSym(n=n,r=2)
problem["obj"]=objF;
problem["grad"]=gradF;
problem["hessian"]=hessF;
problem["retraction"]="proj"
problem["control","tol"]=0.1
problem["control","Delta0"]=0.5
problem["control","DeltaMax"]=10
problem["control","rhoMin"]=0.2
problem["control","iterMax"]=50
problem["control","iterSubMax"]=300
problem["control","conjMethod"]="PR"
lambda=100
problem["control","alpha"]=10
problem["control","beta"]=0.8
problem["control","sigma"]=0.6
res=trustRegion(problem)
res=steepestDescent(problem)
#res=conjugateGradient(problem)


res$NumIter
res$optValue

vHat=res$optY[[1]]
svd(vHat)$d[1:3]
eigen(vHat)$values[1:3]

par(mfrow=c(2,1),mar=c(1,4,0.8,0.5))
pc1=Basis%*%eigen(vHat)$vectors[,1]
plot(x,pc1,type="l",ylab="PC1",lty=2)
lines(x,truFun[,1])
pc2=Basis%*%eigen(vHat)$vectors[,2]
plot(x,pc2,type="l",ylab="PC2",lty=2)
lines(x,-truFun[,2])
legend(0,1.5,c("True","Estimated"),lty=c(1,2))

