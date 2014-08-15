##################################################################################
#This is an eample just to test the special linear group manifold
#In simulation, A is an n*n matrix, x is a special linear transformation of Y
#We then try to find the most similar matrix of A from X

library(expm)
set.seed(120)
n=4
A=matrix(runif(n*n),n,n)

set.seed(88)
B=matrix(runif(n*n),n,n)
B=B-(sum(diag(B))/n)*diag(rep(1,n))
sum(diag(B))
C=expm(B)
det(C)

D=C%*%A

problem=specialLinear(n)
problem["obj"]=function(X){
  sum((A-X%*%D)^2)
}

problem["grad"]=function(X){
  -2*(A-X%*%D)%*%t(D)
}

problem["hessian"]=function(X,Z){
  2*Z%*%D%*%t(D)
}

problem["retraction"]="Norm"
problem["control","tol"]=0.00001
problem["control","Delta0"]=3
problem["control","DeltaMax"]=10
problem["control","rhoMin"]=0.01
problem["control","iterMax"]=1000
problem["control","alpha"]=3
problem["control","iterSubMax"]=100
problem["control","conjMethod"]="PR"
problem["control","threadNum"]=1
problem["control","particleNum"]=200
problem["control","omega"]=0.8
problem["control","phi1"]=0.2
problem["control","phi2"]=0.8

res=steepestDescent(problem)
Y=problem@Y[[1]]
det(Y) #should be 1

H=problem["grad"](Y)
H=H-1/4*sum(diag(H%*%solve(Y)))*Y
sum(diag(H%*%solve(Y))) #should be 0 if correct 
res$optValue   # should be very close to zero if correct
