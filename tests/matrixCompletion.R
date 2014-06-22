set.seed(104)
U=matrix(runif(60),ncol=3)
V=matrix(runif(57),ncol=3)
sigma=diag(c(1,1,1))
lowR=U%*%sigma%*%t(V)

mask=matrix(runif(20*19)>0.3,ncol=19)
sum(mask)


problem=fixedRank(20,19,3)
problem["obj"]=function(X){
  temp=0.5*(lowR-X)^2
  sum(temp[mask])
}
problem["grad"]=function(X){
  temp=-lowR+X
  temp[!mask]=0
  temp
}
problem["hessian"]=function(X,Z){
  temp=Z
  temp[!mask]=0
  temp
}



problem["control","tol"]=0.00001
problem["control","Delta0"]=10
problem["control","DeltaMax"]=20
problem["control","rhoMin"]=0.05
problem["control","iterMax"]=1000
problem["control","iterSubMax"]=1000
problem["control","conjMethod"]="FR"
problem["control","alpha"]=0.05
problem["control","sigma"]=0.01
#res=steepestDescent(problem)
#res=trustRegion(problem)
res=conjugateGradient(problem)
res$optValue
res$NumIter
 temp=res$optY[[1]]-lowR
 sum((temp[!mask])^2)

