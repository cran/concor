\name{svdcp}
\title{SVD for a Column Partitioned matrix x}
\description{SVD for a Column Partitioned matrix x. r global successive solutions}
\keyword{}
\keyword{}
\alias{svdcp}

\usage{
svdcp(x,H,r)
}

\arguments{
\item{x}{is a p x q matrix}
\item{H}{is a row vector which contains the numbers qi, i=1,...,kx, of
  the partition of x with kx column blocks xi : \eqn{\sum q_i = q.}}
\item{r}{is the wanted number of successive solutions.}
}

\value{
list with following components
\item{u}{is a p x r matrix; u'*u = Identity}
\item{v}{is a q x r matrix of kx row blocks vi (qi x r); vi'*vi = Identity}
\item{s2}{is a kx x r matrix; each column k contains kx values
\eqn{(u[,k]'*x_i*v_i[,k])^2}, the partial (squared) singular values
relative to xi}
}

\details{
The first solution calculates 1+kx normed vectors: the vector u[,1] of
\eqn{R^p} associated to the kx vectors vi[,1]'s of \eqn{R^{q_i}}. by maximizing
\eqn{\sum_i (u[,1]'*x_i*v_i[,1])^2}, with 1+kx norm constraints.  A
value \eqn{(u[,1]'*x_i*v_i[,1])^2} measures the relative link between
\eqn{R^p} and \eqn{R^{q_i}} associated to xi. It corresponds to a partial squared
singular value notion, since \eqn{\sum_i (u[,1]'*x_i*v_i[,1])^2=s^2},
where s is the usual first singular value of x.

The second solution is obtained from the same criterion, but after
replacing each xi by xi-xi*vi[,1]*vi[,1]'.  And so on for the
successive solutions 1,2,...,r .  The biggest number of solutions may
be r=inf(p,qi), when the xi's are supposed with full rank; then
rmax=min([min(H),p]).
}

\references{
Lafosse R. & Hanafi M.(1997) Concordance d'un tableau avec K tableaux:
Definition de K+1 uples synthetiques. Revue de Statistique Appliquee
vol.45,n.4.
}

\examples{
x<-matrix(runif(200),10,20)
s<-svdcp(x,c(5,5,10),1)
ss<-svd(x);ss$d[1]^2
sum(s$s2)
}
