#This file provides numerical evaluation of Gradient and Hessian
#They are evaluated from objective function
#They also provide function to check consistency of user-specified functions


setGeneric("checkGradient",function(object){standardGeneric("trustRegion")})
setMethod("checkGradient","manifold",
          definition=function(object){
            if(is.null(object@obj)) steop("The objective function cound not be NULL.")
            if(is.null(object@grad)) steop("The gradient function cound not be NULL.")
            checkGradient2(object)
          })

checkGradient2=function(object){
  n=object@n
  p=object@p
  Y=matrix(10*rnorm(n*p),n,p)
  funGrad=object@grad(Y) #Gradient from function
  numGrad=numericalGradient(object@obj,Y,n,p)
  max(abs(funGrad-numGrad))
}

numericalGradient<-function(objF,Y,n,p){
  epsilon=0.000001
  f0=objF(Y)
  numGrad=matrix(0,n,p) #numerical Gradient
  for(i in 1:n){
    for(j in 1:p){
      Y1=Y
      Y1[i,j]=Y1[i,j]+epsilon
      f1=objF(Y1)
      numGrad[i,j]=(f1-f0)/epsilon
    }
  }
  numGrad
}


setGeneric("checkHessian",function(object){standardGeneric("trustRegion")})
setMethod("checkHessian","manifold",
          definition=function(object){
            if(is.null(object@obj)) steop("The objective function cound not be NULL.")
            if(is.null(object@hessian)) steop("The hessian function cound not be NULL.")
            checkHessian2(object)
          })

checkHessian2=function(object){
  n=object@n
  p=object@p
  Y=matrix(10*rnorm(n*p),n,p)
  Z=matrix(rnorm(n*p),n,p)
  funHess=object@hessian(Y,Z) #Gradient from function
  numHess=numericalHessian(object@obj,Y,Z,n,p)
  max(abs(funHess-numHess))
}

numericalHessian<-function(objF,Y,Z,n,p){
  epsilon=0.01
  numGrad0=numericalGradient(objF,Y,n,p)
  Y=Y+epsilon*Z
  numGrad1=numericalGradient(objF,Y,n,p)
  (numGrad1-numGrad0)/epsilon
}


numericalHessian2<-function(gradF,Y,Z){  # Numerical hessian from user-specified gradient function
  epsilon=0.01
  numGrad0=gradF(Y)
  Y=Y+epsilon*Z
  numGrad1=gradF(Y)
  (numGrad1-numGrad0)/epsilon
}
