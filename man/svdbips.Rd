\name{svdbips}
\title{SVD for bipartitioned matrix x}
\description{SVD for bipartitioned matrix x. SIMULTANEOUS SOLUTIONS. ("simultaneous svdbip")}
\keyword{}
\keyword{}
\alias{svdbips}

\usage{
svdbips(x,K,H,r)
}

\arguments{
\item{x}{is a p x q matrix}
\item{K}{is a row vector which contains the numbers pk, k=1,...,kx, of the
partition of x with kx row blocks : \eqn{\sum_k p_k=p}}
\item{H}{is a row vector which contains the numbers qh, h=1,...,ky, of the
partition of x with ky column blocks : \eqn{\sum_h q_h=q}}
\item{r}{is the wanted number of solutions}
}

\value{
list with following components
\item{u}{is a p x r matrix of kx row blocks uk (pk x r); uk'*uk = Identity}
\item{v}{is a q x r matrix of ky row blocks vh (qh x r); vh'*vh = Identity}
\item{s2}{is a kx x ky x r array; for a fixed solution k, each matrix s2[,,k] contains kxky values \eqn{(u_h'*x_{kh}*v_k)^2}, the "partial (squared) singular values" relative to \eqn{x_{kh}}.}
}

\details{
One set of r solutions is calculated by maximizing \eqn{\sum_i \sum_k \sum_h
(u_k[,i]'*x_{kh}*v_h[,i])^2}, with kx+ky orthonormality constraints (for
each uk and each vh).  For each fixed r value, the solution is totally
new (does'nt consist to complete a previous calculus of one set of r-1
solutions).  rmax=min([min(K),min(H)]).  When r=1, it is svdbip (thus
it is svdcp when r=1 and kx=1).

Convergence of algorithm may be not global. So the below proposed
initialisation of the algorithm may be not very suitable for some data
sets.  Several different random initialisations with normed vectors
might be considered and the best result then choosen....  }

\references{
Lafosse R. & Ten Berge J. A simultaneous CONCOR method for the analysis of two
partitioned matrices. submitted.
}

\examples{
x<-matrix(runif(200),10,20)
s1<-svdbip(x,c(3,4,3),c(5,5,10),2);sum(sum(sum(s1$s2)))
ss<-svdbips(x,c(3,4,3),c(5,5,10),2);sum(sum(sum(ss$s2)))
}

