set.seed(120)
nn=20
A=matrix(runif(nn^2),nn,nn)
A=A+t(A)
problem=grassmannQ(n=nn,p=2)
problem["obj"]=function(X){
  -sum(diag(t(X)%*%A%*%X))
}
#set gradient function
problem["grad"]=function(X){
  -2*A%*%X
}
#set hessian function
problem["hessian"]=function(X,Z){
  -2*A%*%Z
}
problem["retraction"]="QR"
problem["control","tol"]=0.001
problem["control","Delta0"]=5
problem["control","DeltaMax"]=10
problem["control","rhoMin"]=0.01
problem["control","iterMax"]=1000
problem["control","alpha"]=1
problem["control","iterSubMax"]=1000
problem["control","particleNum"]=100
problem["control","conjMethod"]="FR"
checkGradient(problem)
checkHessian(problem)
res=particleSwarm(problem)
-res$optValue
sum(eigen(A)$values[1:2])
res$NumIter

#eigen(A)$vectors[,1:2]
#res$optY[[1]]
#svd(cbind(eigen(A)$vectors[,1:2],res$optY[[1]]))


set.seed(120)
nn=20
A=matrix(runif(nn^2),nn,nn)
A=A+t(A)
problem=grassmannSub(n=nn,r=2)
problem["obj"]=function(X){
  -sum(diag(A%*%X))
}
#set gradient function
problem["grad"]=function(X){
  -A
}
#set hessian function
problem["hessian"]=function(X,Z){
  -matrix(0,nn,nn)
}
problem["retraction"]="Cayley"
problem["control","tol"]=0.1
problem["control","Delta0"]=3
problem["control","DeltaMax"]=10
problem["control","rhoMin"]=0.05
problem["control","iterMax"]=2000
problem["control","alpha"]=5
problem["control","iterSubMax"]=1000
problem["control","particleNum"]=100
problem["control","conjMethod"]="FR"
checkGradient(problem)
checkHessian(problem)
res=conjugateGradient(problem)
-res$optValue
sum(eigen(A)$values[1:2])
res$NumIter


