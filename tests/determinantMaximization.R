######################################################################################
###This is an example to use rOptManifold package to optimize the object function:  ##
###det(B+X'AX), where A is symmetirc, with constraint trace(X'X)=1. We use the class ##
###of sphere manifold. That is equivalent to letting function be -det(B+X'AX) then  ##
###find the minimizer.          ######################################################
###We use several solver to get the optimizer and optimal value.##############
######################################################################################

set.seed(88)
n=20
A=matrix(rnorm(n^2),n,n)
A=t(A)+A
B=matrix(rnorm(2^2),2,2)
problem=sphere(n,p=2)

#enter object function
problem["obj"]=function(X){
  -det(B+t(X)%*%A%*%X)
}

# Minor and cofactor
minor <- function(A, i, j) {det( as.matrix(A[-i,-j]) )}
cofactor <- function(A, i, j) (-1)^(i+j) * minor(A,i,j)

# to get adjugate function with a loop
adjugate <- function(A) {
  n <- nrow(A)
  B <- matrix(NA, n, n)
  for( i in 1:n )
    for( j in 1:n )
      B[j,i] <- cofactor(A, i, j)
  B
}

#enter gradient on ambient space
problem["grad"]=function(X){
  C<-adjugate(B+t(X)%*%A%*%X)
  -(A%*%X%*%t(C)+A%*%X%*%C)
}

#set control parameters and the retraction method of sphere manifold

problem["control","tol"]=0.01
problem["control","Delta0"]=5
problem["control","DeltaMax"]=10
problem["control","rhoMin"]=0.01

problem["control","alpha"]=1
problem["control","iterSubMax"]=100
problem["control","conjMethod"]="FR"

problem["control","threadNum"]=4

problem["control","particleNum"]=100
problem["control","iterMax"]=1000
problem["retraction"]="Exp"
ptm <- proc.time()
res<-steepestDescent(problem)
ptm2 <- proc.time()
#res<-particleSwarm(problem)

res2<-conjugateGradient(problem)
proc.time() - ptm2
ptm2 - ptm

res$optValue
res2$optValue
#sum(res$optY[[1]]^2)
