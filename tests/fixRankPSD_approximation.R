set.seed(120)
nn=40
rr=5
A=matrix(runif(nn*rr),nn,rr)
A=A%*%t(A)
scaling=diag(1/sqrt(diag(A)))
A=scaling%*%A%*%scaling
#scaling=1/sum(diag(A))
#A=A/scaling
B=A#+matrix(rnorm(nn*rr,sd=0.0001),nn,nn)

#problem=fixedRankPSD(n=nn,p=nn,r=rr)
#problem=spectahedron(n=nn,p=nn,r=rr)
problem=elliptope(n=nn,p=nn,r=rr)
problem["obj"]=function(X){
  sum((B-X%*%t(X))^2)
}
#set gradient function
problem["grad"]=function(X){
  4*(X%*%t(X)-B)%*%X
}
#set hessian function
problem["hessian"]=function(X,Z){
  4*(X%*%t(Z)+Z%*%t(X))%*%X+4*(X%*%t(X)-B)%*%Z
}
problem["retraction"]="proj"
problem["control","tol"]=0.001
problem["control","Delta0"]=3
problem["control","DeltaMax"]=8
problem["control","rhoMin"]=0.1
problem["control","iterMax"]=1000
problem["control","alpha"]=5
problem["control","iterSubMax"]=1000
problem["control","conjMethod"]=1

res=trustRegion(problem)
res=steepestDescent(problem)
res=conjugateGradient(problem)
res$optValue
aHat=(res$optY[[1]])
aHat=aHat%*%t(aHat)
sum(abs(aHat-A))
aHat[1:3,1:3]
A[1:3,1:3]

res$NumIter


