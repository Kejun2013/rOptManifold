set.seed(110)
pp=100#dimension p
nn=50#number of obs
rr=1
pca1=matrix(c(1:5,rep(0,pp-5)),ncol=1)
pca1=pca1/norm(pca1,"F")
Y=pca1%*%matrix(rnorm(nn,0,5),nrow=1)+matrix(rnorm(pp*nn,0,0.5),nrow=pp)
Sigma=cov(t(Y))

##tuning parameter
rho=0.2

problem=spectahedron(n=pp,r=rr)
#problem=elliptope(n=nn,p=nn,r=rr)
problem["obj"]=function(X){
  XXt=X%*%t(X)
  sparsity=sum(sqrt(X^2+0.01))
  -sum(diag(Sigma%*%XXt))+rho*sparsity
}
#set gradient function
problem["grad"]=function(X){
  sparsity=X/(sqrt(X^2+0.005))
  gradF=-2*Sigma%*%X+rho*sparsity
  gradF
}
#set hessian function
problem["hessian"]=function(X,Z){
  sparsity1=Z/(sqrt(X^2+0.005))
  sparsity2=X^2*sparsity1/(X^2+0.005)
  eucH=-2*Sigma%*%Z+rho*sparsity1-rho*sparsity2
  eucH
}
problem["retraction"]="proj"
problem["control","tol"]=0.01
problem["control","Delta0"]=10
problem["control","DeltaMax"]=10
problem["control","rhoMin"]=0.2
problem["control","iterMax"]=500
problem["control","alpha"]=2
problem["control","iterSubMax"]=1000
problem["control","conjMethod"]="PR"
problem["control","threadNum"]=1
problem["control","particleNum"]=200
problem["control","omega"]=0.8
problem["control","phi1"]=1
problem["control","phi2"]=1.5

#res=trustRegion(problem)
res=steepestDescent(problem)
#res=conjugateGradient(problem)
#res=particleSwarm(problem)
res$NumIter
res$optValue
pcaHat=res$optY[[1]]
sum(pcaHat^2)
pcaHat[abs(pcaHat)<0.05]=0
sum(abs(pcaHat)>0)
pcaHat[1:10]
pca1[1:10]


