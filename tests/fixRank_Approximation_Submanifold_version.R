set.seed(120)
nn=40
rr=3
A=matrix(runif(nn*rr),nn,rr)
A=A%*%t(A)

#add noise
B=A#+matrix(rnorm(nn*rr,sd=0.05),nn,nn)

problem=fixedRankSym(n=nn,r=rr)
problem["obj"]=function(X){
  sum((B-X)^2)
}
#set gradient function
problem["grad"]=function(X){
  2*(X-B)
}
#set hessian function
problem["hessian"]=function(X,Z){
  2*Z
}
problem["retraction"]="proj"
problem["control","tol"]=0.01
problem["control","Delta0"]=0.5
problem["control","DeltaMax"]=10
problem["control","rhoMin"]=0.05
problem["control","iterMax"]=1000
problem["control","alpha"]=10
problem["control","iterSubMax"]=100
problem["control","conjMethod"]="PR"

res=trustRegion(problem)
res=steepestDescent(problem)
res=conjugateGradient(problem)
res$optValue
res$NumIter

aHat=(res$optY[[1]])
max(abs(aHat-A))
aHat[1:3,1:3]
A[1:3,1:3]


