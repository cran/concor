\name{concorcano}
\title{Canonical analysis of several sets with another set}
\description{Relative proximities of several subsets of variables Yj with another set X. SUCCESSIVE SOLUTIONS}
\keyword{}
\keyword{}
\alias{concorcano}

\usage{
concorcano(x,y,py,r)
}

\arguments{
\item{x}{is a n x p matrix of p centered variables}
\item{y}{is a n x q matrix of q centered variables}
\item{py}{is a row vector which contains the numbers qi, i=1,...,ky, of
  the ky subsets yi of y : \eqn{\sum_i q_i}=sum(py)=q. py is the partition vector of y}
\item{r}{is the wanted number of successive solutions}
}

\value{
list with following components
\item{cx}{is n x r matrix of the r canonical components of x}
\item{cy}{is n.ky x r matrix. The ky blocks cyi of the rows n*(i-1)+1 : n*i contain the r canonical components relative to Yi}
\item{rho2}{is a ky x r matrix; each column k contains ky squared
  canonical correlations \eqn{\rho(cx[,k],cy_i[,k])^2}}
}

\details{
The first solution calculates a standardized canonical component
cx[,1] of x associated to ky standardized components cyi[,1] of yi by
maximizing \eqn{\sum_i \rho(cx[,1],cy_i[,1])^2}.

The second solution is obtained from the same criterion, with ky
orthogonality constraints for having rho(cyi[,1],cyi[,2])=0 (that
implies rho(cx[,1],cx[,2])=0).  For each of the 1+ky sets, the r
canonical components are 2 by 2 zero correlated.

The ky matrices (cx)'*cyi are triangular.

This function uses concor function.
}

\references{
Hanafi & Lafosse (2001) Generalisation de la regression lineaire
simple pour analyser la dependance de K ensembles de variables avec un
K+1 eme.  Revue de Statistique Appliquee vol.49, n.1
}

\examples{
x<-matrix(runif(50),10,5);y<-matrix(runif(90),10,9)
x<-scale(x);y<-scale(y)
ca<-concorcano(x,y,c(3,2,4),2)
diag(t(ca$cx)\%*\%ca$cy[1:10,]/10)^2
ca$rho2[1,]
}

