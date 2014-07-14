#This file provides numerical evaluation of Gradient and Hessian
#They are evaluated from objective function
#They also provide function to check consistency of user-specified functions


checkGradient2=function(object){
  n=object@n
  p=object@p
  if(length(n)==1){
    Y=matrix(10*rnorm(n*p),n,p)
    funGrad=object@grad(Y) #Gradient from function
    numGrad=numericalGradient(object@obj,Y,n,p,0)
    return(max(abs(funGrad-numGrad)))
  }else{
    numError=rep(0,length(n))
    Y=list()
    for(k in 1:length(n)){
      Y=c(Y,list(matrix(10*rnorm(n[k]*p[k]),n[k],p[k])))
    }
    for(k in 1:length(n)){
      funGrad=object@grad(Y,k) #Gradient from function
      numGrad=numericalGradient(object@obj,Y,n,p,k)
      numError[k]=max(abs(funGrad-numGrad))
    }
    return(max(numError))
  }
}

numericalGradient<-function(objF,Y,n,p,k){
  epsilon=0.0001
  f0=objF(Y)
  if(k==0){
    numGrad=matrix(0,n,p) #numerical Gradient
    for(i in 1:n){
      for(j in 1:p){
        Y1=Y
        Y1[i,j]=Y1[i,j]+epsilon
        f1=objF(Y1)
        numGrad[i,j]=(f1-f0)/epsilon
      }
    }
  }else{
    numGrad=matrix(0,n[k],p[k]) #numerical Gradient
    for(i in 1:n[k]){
      for(j in 1:p[k]){
        Y2=Y
        Y1=Y[[k]]
        Y1[i,j]=Y1[i,j]+epsilon
        Y2[[k]]=Y1
        f1=objF(Y2)
        numGrad[i,j]=(f1-f0)/epsilon
      }
    }
  }
  numGrad
}


checkHessian2=function(object){
  n=object@n
  p=object@p
  if(length(n)==1){
    Y=matrix(10*rnorm(n*p),n,p)
    Z=matrix(rnorm(n*p),n,p)
    funHess=object@hessian(Y,Z) #Gradient from function
    numHess=numericalHessian(object@obj,Y,Z,n,p)
    return(max(abs(funHess-numHess)))
  }else{
    numError=rep(0,length(n))
    Y=list()
    for(k in 1:length(n)){
      Y=c(Y,list(matrix(10*rnorm(n[k]*p[k]),n[k],p[k])))
    }
    for(k in 1:length(n)){
      Z=matrix(rnorm(n[k]*p[k]),n[k],p[k])
      funHess=object@hessian(Y,Z,k) #Gradient from function
      numHess=numericalHessian(object@obj,Y,Z,n,p,k)
      numError[k]=max(abs(funHess-numHess))
    }
    return(max(numError))
  }
}

numericalHessian<-function(objF,Y,Z,n,p,k=0){
  epsilon=0.001
  if(k==0){
    numGrad0=numericalGradient(objF,Y,n,p,k)
    Y=Y+epsilon*Z
    numGrad1=numericalGradient(objF,Y,n,p,k)
    return((numGrad1-numGrad0)/epsilon)
  }else{
    numGrad0=numericalGradient(objF,Y,n,p,k)
    Y[[k]]=Y[[k]]+epsilon*Z
    numGrad1=numericalGradient(objF,Y,n,p,k)
    return((numGrad1-numGrad0)/epsilon)
  }
}


numericalHessian2<-function(gradF,Y,Z,k=0){  # Numerical hessian from user-specified gradient function
  epsilon=0.001
  if(k>0){
    numGrad0=gradF(Y,k)
    Y[[k]]=Y[[k]]+epsilon*Z
    numGrad1=gradF(Y,k)
    return((numGrad1-numGrad0)/epsilon)
  }else{
    numGrad0=gradF(Y)
    Y=Y+epsilon*Z
    numGrad1=gradF(Y)
    return((numGrad1-numGrad0)/epsilon)
  }
}
