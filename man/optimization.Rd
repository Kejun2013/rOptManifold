\name{optimization}
\alias{steepestDescent}
\alias{trustRegion}
\alias{conjugateGradient}
\alias{particleSwarm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimization Algorithms
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
steepestDescent(object)

trustRegion(object)

conjugateGraident(object)

particleSwarm(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of manifold class. See \link{stiefel}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{optY }{Optimizer at which function is maximized}
\item{NumIter }{Number of iteration used until convergence}
\item{optValue}{The minimum value achieved by the objective function}
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
##A graph cut problem
n=40##number of vertices
set.seed(50)
B=C=D=matrix(0,n/2,n/2)
##The graph is devided into two group
##First group is 1:(n/2), the second groupd is (n/2+1):n
##Within each group, with is 0.2 prob that an edge exsits
##Between group, the prob is 0.6
fHalf=1:(n/2);sHalf=(n/2+1):n
elSel=upper.tri(B)
B[elSel]=(runif(sum(elSel))<0.03)
C[elSel]=(runif(sum(elSel))<0.3)
D[elSel]=(runif(sum(elSel))<0.05)
 
# ##A is adjacency matrix and L is Laplacian matrix
 A=matrix(0,n,n)
 A[fHalf,fHalf]=B;A[fHalf,sHalf]=C;A[sHalf,sHalf]=D;
 A=A+t(A)
 w=colSums(A)
 L=diag(w)-A
 
##plot the graph
 par(mar=rep(0,4))
 plot(c(1,5),c(1,(n/4+1)),type="n",bty="n",xaxt="n",yaxt="n")
 px=c(rep(1.6,n/4),rep(2,n/4),rep(4,n/4),rep(4.6,n/4))
 py=c(1:(n/4)+0.5,1:(n/4),1:(n/4)+0.5,1:(n/4))
 
 for(i in 1:n){
   for(j in 1:n){
     if(A[i,j]==1){
       segments(px[i],py[i],px[j],py[j])
     }
   }
 }
 
 for(i in 1:n){
   points(px[i],py[i],pch=16,col="grey",cex=3)
}
}

