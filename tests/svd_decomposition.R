library(rOptManifold)
#problem_svd=stiefel(n=3,p=2)*stiefel(n=4,p=2)
problem_svd=stiefel(n=c(3,4),p=c(2,2))


B=matrix(c(3,1,2,1,3,1,1,1,1,1,1,1),nrow=3);B
N=diag(c(2,1));N
problem_svd["obj"]=function(X){
   -sum(diag(t(X[[1]])%*%B%*%X[[2]]%*%N))
}


problem_svd["grad"]=function(X,k){
  if(k==1){
    return(-B%*%X[[2]]%*%N)
  }else{
     return(-t(B)%*%X[[1]]%*%N)
  }
}

problem_svd["hessian"]=function(X,Z,k){
  if(k==1){
    n1=dim(X[[1]])[1];p1=dim(X[[1]])[2]
    return(matrix(0,n1,p1))
  }else{
    n1=dim(X[[2]])[1];p1=dim(X[[2]])[2]
    return(matrix(0,n1,p1))
  }
}

problem_svd["retraction"]=c("cayley","cayley") #Exp QR Cayley
checkGradient(problem_svd)
checkHessian(problem_svd)
res=trustRegion(problem_svd)

