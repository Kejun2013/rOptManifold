set.seed(120)
nn=10
A=matrix(runif(nn^2),nn,nn)
A=A+t(A)
N=diag(c(10,5))
problem=stiefel(n=nn,p=2)
problem["obj"]=function(X){
  -sum(diag(t(X)%*%A%*%X%*%N))
}
#set gradient function
problem["grad"]=function(X){
  -2*A%*%X%*%N
}
#set hessian function
problem["hessian"]=function(X,Z){
  -2*A%*%Z%*%N
}
problem["retraction"]="QR"
problem["control","tol"]=0.01
problem["control","Delta0"]=5
problem["control","DeltaMax"]=10
problem["control","rhoMin"]=0.01
problem["control","iterMax"]=1000
problem["control","alpha"]=1
problem["control","iterSubMax"]=1000
problem["control","conjMethod"]="FR"

res=particleSwarm(problem)
-res$optValue
res$NumIter
eigen(A)$values[1:2]%*%diag(N)
