\name{concorgm}
\title{Analyzing a set of partial links between Xi and Yj}
\description{Analyzing a set of partial links between Xi and Yj, SUCCESSIVE SOLUTIONS}
\keyword{}
\keyword{}
\alias{concorgm}

\usage{
concorgm(x,px,y,py,r)
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
\item{u}{is a p x r matrix of kx row blocks ui (pi x r), the orthonormed partial axes of xi; associated partial components: xi*ui}
\item{v}{is a q x r matrix of ky row blocks vj (qj x r), the orthonormed partial axes of yj; associated partial components: yj*vj}
\item{cov2}{is a kx x ky x r array; for r fixed to k, the matrix contains kxky  squared covariances \eqn{\mbox{cov2}(x_i*u_i[,k],y_j*v_j[,k])^2}, the partial links between xi and yj measured with the solution k.}
}

\details{
For the first solution, \eqn{\sum_i \sum_j \mbox{cov2}(x_i*u_i[,1],y_j*v_j[,1])} is the
optimized criterion. The second solution is calculated from the same
criterion, but with \eqn{x_i-x_i*u_i[,1]*u_i[,1]'} and \eqn{y_j-y_j*v_j[,1]*v_j[,1]'}
instead of the kx+ky matrices xi and yj. And so on for the other
solutions.  When kx=1 (px=p), take concor.m

This function uses the svdbip function.
}

\references{
Kissita, Cazes, Hanafi & Lafosse (2004) Deux methodes d'analyse
factorielle du lien entre deux tableaux de variables
partitionnées. Revue de Statistique Appliquée, Vol 52, n° 3, 73-92.
}

\examples{
x<-matrix(runif(50),10,5);y<-matrix(runif(90),10,9)
x<-scale(x);y<-scale(y)
cg<-concorgm(x,c(2,3),y,c(3,2,4),2)
diag(t(x[,1:2]\%*\%cg$u[1:2,])\%*\%y[,1:3]\%*\%cg$v[1:3,]/10)^2
cg$cov2[1,1,]
}

