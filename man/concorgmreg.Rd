\name{concorgmreg}
\title{Regression of subsets Yj by subsets Xi}
\description{Regression of subsets Yj by subsets Xi for comparing all the explanatory-explained pairs (Xi,Yj). SUCCESSIVE SOLUTIONS}
\keyword{}
\keyword{}
\alias{concorgmreg}

\usage{
concorgmreg(x,px,y,py,r)
}

\arguments{
\item{x}{is a n x p matrix of p centered variables}
\item{y}{is a n x q matrix of q centered variables}
\item{px}{is a row vector which contains the numbers pi, i=1,...,kx, of
  the kx subsets xi of x : \eqn{\sum p_i}=sum(px)=p. px is the partition
  vector of the columns of x.}
\item{py}{is the partition vector of y with ky subsets yj, j=1,...,ky. sum(py)=q}
\item{r}{is the wanted number of successive solutions rmax <= min(min(px),min(py),n)}
}

\value{
list with following components
\item{cx}{is a n.kx x r matrix of kx row blocks cxi (n x r). Each row block contains r partial explanatory components}
\item{v}{is a q x r matrix of ky row blocks vj (qj x r), the orthonormed partial axes of yj; The components yj*vj are the explained components}
\item{varexp}{is a kx x ky x r array; for a fixed solution k, the matrix varexp[,,k] contains kxky explained variances obtained by a simultaneous regression of all the yj by all the xi, so the values \eqn{\mbox{rho2}(cx[n*(i-1)+1:n*i,k],y_j*v_j[,k]) var(y_j*v_j[,k])}}
}

\details{
For the first solution, \eqn{\sum_i \sum_j \mbox{rho2}(cx_i[,1],y_j*v_j[,1])
\mbox{var}(y_j*v_j[,1])} is the optimized criterion. The second solution is
calculated from the same criterion, but with \eqn{y_j-y_j*v_j[,1]*v_j[,1]'}
instead of the matrices yj and with orthogonalities for having two by
two zero correlated the explanatory components defined for each matrix
xi. And so on for the other solutions. One solution k associates kx
explanatory components (in cx[,k]) to ky explained components. When
kx =1 (px=p), take concoreg function

This function uses the concorgm function
}

\references{
Hanafi & Lafosse (2004) Regression of a multi-set by another based on an
extension of the SVD. COMPSTAT'2004 Symposium
}

\examples{
x<-matrix(runif(50),10,5);y<-matrix(runif(90),10,9)
x<-scale(x);y<-scale(y)
cr<-concorgmreg(x,c(2,3),y,c(3,2,4),2)
diag(t(cr$cx[1:10,])\%*\%y[,1:3]\%*\%cr$v[1:3,]/10)^2
cr$varexp[1,1,]
}

