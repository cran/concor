\name{concor}
\title{Relative links of several subsets of variables}
\description{Relative links of several subsets of variables Yj with another set X. SUCCESSIVE SOLUTIONS}
\keyword{}
\keyword{}
\alias{concor}

\usage{
concor(x,y,py,r)
}

\arguments{
\item{x,y}{are n x p and n x q matrices of p and q centered columns}
\item{py}{is a row vector which contains the numbers qi, i=1,...,ky, of the ky subsets yi of y : sum(qi)=sum(py)=q. py is the partition vector of y}
\item{r}{is the wanted number of successive solutions}
}

\value{
list with following components
\item{u}{is a p x r matrix of axes in Rp relative to x; u'*u = Identity}
\item{v}{is a q x r matrix of ky row blocks vi (qi x r) of axes in Rqi relative to yi; vi'*vi = Identity}
\item{V}{is a q x r matrix of axes in Rq relative to y; V'*V = Identity}
\item{cov2}{is a ky x r matrix; each column k contains ky squared covariances \eqn{\mbox{cov}(x*u[,k],y_i*v_i[,k])^2}, the partial measures of link}
}

\details{
The first solution calculates 1+kx normed vectors: the vector u[:,1]
of Rp associated to the ky vectors vi[:,1]'s of Rqi, by maximizing
\eqn{\sum_i \mbox{cov}(x*u[,k],y_i*v_i[,k])^2}, with 1+ky norm constraints on the
axes. A component x*u[,k] is associated to ky partial components
yi*vi[,k] and to a global component y*V[,k].  \eqn{\mbox{cov}(x*u[,k],y*V[,k])^2 =
\sum \mbox{cov}(x*u[,k],y_i*v_i[,k])^2}. y*V[,k] is a global component of the
components yi*vi[,k].

The second solution is obtained from the same criterion, but after
replacing each yi by \eqn{y_i-y_i*v_i[,1]*v_i[,1]'}.  And so on for the
successive solutions 1,2,...,r.  The biggest number of solutions may
be r=inf(n,p,qi), when the x'*yi's are supposed with full rank; then
rmax=min(c(min(py),n,p)).  For a set of r solutions, the matrix u'X'YV
is diagonal and the matrices u'X'Yjvj are triangular (good partition
of the link by the solutions).  concor.m is the svdcp.m function
applied to the matrix x'y.
}

\references{
Lafosse R. & Hanafi M.(1997) Concordance d'un tableau avec K tableaux:
Definition de K+1 uples synthetiques. Revue de Statistique Appliquee
vol.45,n.4.
}

\examples{
# To make some "GPA" : so, by posing the compromise X = Y,
# "procrustes" rotations to the "compromise X" then are :
# Yj*(vj*u').

x<-matrix(runif(50),10,5);y<-matrix(runif(90),10,9)
x<-scale(x);y<-scale(y)
co<-concor(x,y,c(3,2,4),2)
((t(x\%*\%co$u[,1])\%*\%y[,1:3]\%*\%co$v[1:3,1])/10)^2;co$cov2[1,1] 
t(x\%*\%co$u)\%*\%y\%*\%co$V
}
