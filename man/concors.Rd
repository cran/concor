\name{concors}
\title{"simultaneous  concorgm"}
\description{concorgm with the set of r solutions simultaneously optimized}
\keyword{}
\keyword{}
\alias{concors}

\usage{
concors(x,px,y,py,r)
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
\item{u}{is a p x r matrix of kx row blocks ui (pi x r), the orthonormed partial axes of xi; associated partial components: xi*ui}
\item{v}{is a q x r matrix of ky row blocks vj (qj x r), the orthonormed partial axes of yj; associated partial components: yj*vj}
\item{cov2}{is a kx x ky x r array; for r fixed to k, the matrix contains kxky  squared covariances \eqn{\mbox{cov}(x_i*u_i[,k],y_j*v_j[,k])^2}, the partial links between xi and yj measured with the solution k}
}

\details{
This function uses the svdbips function
}

\references{See svdbips}

\examples{
x<-matrix(runif(50),10,5);y<-matrix(runif(90),10,9)
x<-scale(x);y<-scale(y)
cs<-concors(x,c(2,3),y,c(3,2,4),2)
diag(t(x[,1:2]\%*\%cs$u[1:2,])\%*\%y[,1:3]\%*\%cs$v[1:3,]/10)^2
cs$cov2[1,1,]
}

