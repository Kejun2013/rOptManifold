set.seed(104)
m = 60
n = 57
U=matrix(runif(m*3),ncol=3)
V=matrix(runif(n*3),ncol=3)
sigma=diag(c(1,1,1))
lowR=U%*%sigma%*%t(V)

mask=matrix(runif(m*n)>0.3,ncol=n)
sum(mask)


problem=fixedRank(m,n,3)
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



problem["control","tol"]=0.001
problem["control","Delta0"]=10
problem["control","DeltaMax"]=20
problem["control","rhoMin"]=0.05
problem["control","iterMax"]=1500
problem["control","iterSubMax"]=1000
problem["control","conjMethod"]="FR"
problem["control","alpha"]=1
problem["control","sigma"]=0.01
problem["control","threadNum"]=1
problem["control","particleNum"]=900
problem["control","omega"]=0.8
problem["control","phi1"]=0.2
problem["control","phi2"]=0.8
res=steepestDescent(problem)
#res=trustRegion(problem)

#res=conjugateGradient(problem)

#res=particleSwarm(problem)
res$optValue
res$NumIter
 temp=res$optY[[1]]-lowR
 max((temp[!mask])^2)

