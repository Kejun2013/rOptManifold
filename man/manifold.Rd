\name{manifold}
\alias{stiefel}
\alias{grassmannSub}
\alias{grassmannQ}
\alias{oblique}
\alias{elliptope}
\alias{fixedRank}
\alias{fixedRankSym}
\alias{fixedRankPSD}
\alias{spectahedron}
\alias{sphere}
\alias{projective}
\alias{specialLinear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Manifold
}
\description{
Functions for creating manifold of various type. 
}
\usage{
sphere(n,p,retraction="Norm")

stiefel(n,p,retraction="QR")

grassmannQ(n,p,retraction="QR")

grassmannSub(n,r,retraction="QR")

fixedRank(n,p,r, retraction="Proj")

fixedRankSym(n,r,retraction="Proj")

fixedRankPSD(n,r,retraction="Proj")

spectahedron(n,r,retraction="Proj")

elliptope(n,r,retraction="Proj")

oblique(n,p,retraction="Norm")

specialLinear(n,retraction="Norm")

projective(n,p,retraction="Norm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
   A scalar or vector containing row number(s) of matrix. When n is a vector of length greater than one, a product manifold is created.
}
  \item{p}{A scalar or vector of the same length as n. The number(s) of columns of matrix(ces).}
  \item{r}{A scalar or vector of the same length as n. The rank(s) of matrix(ces).}
  \item{retraction}{Retraction methods to be used in optimization. See details.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# a simple example of sparse PCA
set.seed(110)
pp=100#dimension p
nn=50#number of obs
rr=1
pca1=matrix(c(1:5,rep(0,pp-5)),ncol=1)
pca1=pca1/norm(pca1,"F")
Y=pca1%*%matrix(rnorm(nn,0,5),nrow=1)+matrix(rnorm(pp*nn,0,0.5),nrow=pp)
Sigma=cov(t(Y))
 
##tuning parameter
rho=0.1
 
problem=spectahedron(n=pp,r=rr)
#problem=elliptope(n=nn,p=nn,r=rr)
problem["obj"]=function(X){
  XXt=X%*%t(X)
  sparsity=sum(sqrt(X^2+0.01))
  -sum(diag(Sigma%*%XXt))+rho*sparsity
}
#set gradient function
problem["grad"]=function(X){
  sparsity=X/(sqrt(X^2+0.01))
  gradF=-2*Sigma%*%X+rho*sparsity
  gradF
}
#set hessian function
problem["hessian"]=function(X,Z){
  sparsity1=Z/(sqrt(X^2+0.001))
  sparsity2=X^2*sparsity1/(X^2+0.001)
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
 
res=trustRegion(problem)
#res=steepestDescent(problem)
#res=conjugateGradient(problem)
#res=particleSwarm(problem)
res$NumIter
res$optValue
pcaHat=res$optY[[1]]
sum(pcaHat^2)
pcaHat[abs(pcaHat)<0.05]=0
sum(pcaHat>0)
pcaHat[1:10]
pca1[1:10]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
