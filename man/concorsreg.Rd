\name{concorsreg}
\title{"simultaneous concorgmreg"}
\description{concorgmreg with the set of r solutions simultaneously optimized}
\keyword{}
\keyword{}
\alias{concorsreg}

\usage{
concorsreg(x,px,y,py,r)
}

\arguments{
\item{x}{is a n x p matrix of p centered variables}
\item{y}{is a n x q matrix of q centered variables}
\item{px}{is a row vector which contains the numbers pi, i=1,...,kx, of the kx subsets xi of x : sum(pi)=sum(px)=p. px is the partition vector of x}
\item{py}{is the partition vector of y with ky subsets yj, j=1,...,ky}
\item{r}{is the wanted number of successive solutions rmax <= min(min(px),min(py),n)}
}

\value{
list with following components
\item{cx}{is a n.kx x r matrix of kx row blocks cxi (n x r). Each row block contains r partial explanatory components}
\item{v}{is a q x r matrix of ky row blocks vj (qj x r), the orthonormed partial axes of yj; The components yj*vj are the explained components.}
\item{varexp}{is a kx x ky x r array; for a fixed solution k, the matrix
 varexp[,,k] contains kxky explained variances obtained by a simultaneous
 regression of all the yj by all the xi, so the values \eqn{\mbox{rho2}(cx[n*(i-1)+1:n*i,k],y_j*v_j[,k]) var(y_j*v_j[,k])}}
}

\details{
This function uses the concors function
}

\references{See svdbips}

\examples{
x<-matrix(runif(50),10,5);y<-matrix(runif(90),10,9)
x<-scale(x);y<-scale(y)
crs<-concorsreg(x,c(2,3),y,c(3,2,4),2)
diag(t(crs$cx[1:10,])\%*\%y[,1:3]\%*\%crs$v[1:3,]/10)^2
crs$varexp[1,1,]
}

