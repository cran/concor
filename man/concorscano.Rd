\name{concorscano}
\title{"simultaneous concorgmcano"}
\description{concorgmcano with the set of r solutions simultaneously optimized}
\keyword{}
\keyword{}
\alias{concorscano}

\usage{
concorscano(x,px,y,py,r)
}

\arguments{
\item{x}{is a n x p matrix of p centered variables}
\item{y}{is a n x q matrix of q centered variables}
\item{px}{is a row vector which contains the numbers pi, i=1,...,kx, of
  the kx subsets xi of x : \eqn{\sum_i p_i}=sum(px)=p. px is the partition vector of x}
\item{py}{is the partition vector of y with ky subsets yj, j=1,...,ky}
\item{r}{is the wanted number of successive solutions rmax <= min(min(px),min(py),n)}
}


\value{
list with following components
\item{cx}{is a n.kx x r matrix of kx row blocks cxi (n x r). Each row block contains r partial canonical components} 
\item{cy}{is a n.ky x r matrix of ky row blocks cyj (n x r). Each row block contains r partial canonical components}
\item{rho2}{is a kx x ky x r array; for a fixed solution k, rho2[,,k] contains kxky squared correlations \eqn{\rho(cx[n*(i-1)+1:n*i,k],cy[n*(j-1)+1:n*j,k])^2}, simultaneously calculated between all the yj with all the xi}
}

\details{
This function uses the concors function
}

\references{See svdbips}

\examples{
x<-matrix(runif(50),10,5);y<-matrix(runif(90),10,9)
x<-scale(x);y<-scale(y)
cca<-concorscano(x,c(2,3),y,c(3,2,4),2)
diag(t(cca$cx[1:10,])\%*\%cca$cy[1:10,]/10)^2
cca$rho2[1,1,]
}

