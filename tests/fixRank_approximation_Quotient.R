set.seed(120)
nn=40
rr=3
A=matrix(runif(nn*rr),nn,rr)
A=A%*%t(A)

#add noise
B=A+matrix(rnorm(nn*rr,sd=0.05),nn,nn)

problem=fixedRankPSD(n=nn,r=rr)
#problem=spectahedron(n=nn,p=nn,r=rr)
#problem=elliptope(n=nn,p=nn,r=rr)
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
problem["control","tol"]=0.1
problem["control","Delta0"]=1
problem["control","DeltaMax"]=5
problem["control","rhoMin"]=0.01
problem["control","iterMax"]=1000
problem["control","alpha"]=1
problem["control","iterSubMax"]=2000
problem["control","conjMethod"]="PR"
problem["control","threadNum"]=1
problem["control","particleNum"]=800
problem["control","omega"]=0.9
problem["control","phi1"]=2
problem["control","phi2"]=2

#res=trustRegion(problem)
res=steepestDescent(problem)
#res=conjugateGradient(problem)
#res=particleSwarm(problem)
res$optValue
res$NumIter

aHat=(res$optY[[1]])
aHat=aHat%*%t(aHat)
max(abs(aHat-A))
aHat[1:3,1:3]
A[1:3,1:3]



