set.seed(120)
nn=40
rr=3
A=matrix(runif(nn*rr),nn,rr)
A=A%*%t(A)

#add noise
B=A+matrix(rnorm(nn*rr,sd=0.05),nn,nn)

problem=projective(n=nn,p=rr)
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
problem["retraction"]="Norm"
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
#res=steepestDescent(problem)
res=conjugateGradient(problem)
#res=particleSwarm(problem)
res$optValue
res$NumIter

aHat=(res$optY[[1]])
aHat=aHat%*%t(aHat)
max(abs(aHat-A))
aHat[1:3,1:3]
A[1:3,1:3]




#############################
#compare
set.seed(100)
nn=50
rr=4
A=matrix(runif(nn*rr),nn,rr)
A=A%*%t(A)

#add noise
B=A+matrix(rnorm(nn*rr,sd=0.05),nn,nn)

problem=projective(n=nn,p=rr)
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
problem["retraction"]="Norm"
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

steepestMethod_pro=function(){
  problem@Y[[1]]=diag(1,50,4)
  res=steepestDescent(problem)
}
trustRegionMethod_pro=function(){
  problem@Y[[1]]=diag(1,50,4)
  res=trustRegion(problem)
}


conjugateMethod_pro=function(){
  problem@Y[[1]]=diag(1,50,4)
  res=conjugateGradient(problem)
}

problem2=fixedRankPSD(n=nn,r=rr)
problem2["obj"]=function(X){
  sum((B-X%*%t(X))^2)
}
#set gradient function
problem2["grad"]=function(X){
  4*(X%*%t(X)-B)%*%X
}
#set hessian function
problem2["hessian"]=function(X,Z){
  4*(X%*%t(Z)+Z%*%t(X))%*%X+4*(X%*%t(X)-B)%*%Z
}
problem2["retraction"]="proj"
problem2["control","tol"]=0.1
problem2["control","Delta0"]=1
problem2["control","DeltaMax"]=5
problem2["control","rhoMin"]=0.01
problem2["control","iterMax"]=1000
problem2["control","alpha"]=1
problem2["control","iterSubMax"]=2000
problem2["control","conjMethod"]="PR"
problem2["control","threadNum"]=1
problem2["control","particleNum"]=800
problem2["control","omega"]=0.9
problem2["control","phi1"]=2
problem2["control","phi2"]=2

steepestMethod_fix=function(){
  problem2@Y[[1]]=diag(1,50,4)
  res=steepestDescent(problem2)
}
trustRegionMethod_fix=function(){
  problem2@Y[[1]]=diag(1,50,4)
  res=trustRegion(problem2)
}


conjugateMethod_fix=function(){
  problem2@Y[[1]]=diag(1,50,4)
  res=conjugateGradient(problem2)
}

require("rbenchmark")
benchmark(steepestMethod_pro(),
          trustRegionMethod_pro(),
          steepestMethod_fix(),
          trustRegionMethod_fix(),
          conjugateMethod_pro(),
          conjugateMethod_fix(),
          order="relative", replications=20)[,1:4]
