
######################################################################################
###This is an example to use rOptManifold package to optimize the procrustes problem:  
###-||AX-B||^2, with contraint diag(X'X)=I_p, i.e. the class of oblique manifold.         
###We use steepestDescent solver to get the optimizer and optimal value, and then 
###test that the gradient of submanifold at the local minimal point is truly zero.
######################################################################################

set.seed(64)
l=10
p=20
n=30
A=matrix(runif(l*n),l,n)
B=matrix(runif(l*p),l,p)
problem=oblique(n,p)

problem["obj"]=function(X){
  C=A%*%X-B
  -sum(diag(t(C)%*%C))
}

#set gradient function
problem["grad"]=function(X){
  -2*t(A)%*%(A%*%X-B)
}

#set hessian function
problem["hessian"]=function(X,Z){
  -2*t(A)%*%A%*%Z
}

problem["retraction"]="NORM"
problem["control","tol"]=0.01
problem["control","Delta0"]=5
problem["control","DeltaMax"]=10
problem["control","rhoMin"]=0.01
problem["control","iterMax"]=1000
problem["control","alpha"]=1
problem["control","iterSubMax"]=1000
problem["control","particleNum"]=100
problem["control","conjMethod"]="FR"
res=steepestDescent(problem)

M=res$optY[[1]]
grad_general=problem["grad"](M)
temp=diag(1,20,20)
diag(temp)=diag(t(M)%*%grad_general)
temp=M%*%temp

#grad_manifold is the gradient of manifold at point M=res$optY[[1]]
grad_manifold=grad_general-temp
grad_manifold
