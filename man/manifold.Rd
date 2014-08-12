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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
