\name{svdbip2}
\title{SVD for bipartitioned matrix x}
\description{SVD for bipartitioned matrix x. r successive Solutions. As SVDBIP, but with another algorithm and another initialisation}
\keyword{}
\keyword{}
\alias{svdbip2}

\usage{
svdbip2(x,K,H,r)
}

\arguments{
\item{x}{is a p x q matrix}
\item{K}{is a row vector which contains the numbers pk, k=1,...,kx, of
  the partition of x with kx row blocks : \eqn{\sum_k p_k=p}}
\item{H}{is a row vector which contains the numbers qh, h=1,...,ky, of the
partition of x with ky column blocks : \eqn{\sum q_h=q}}
\item{r}{is the wanted number of successive solutions}
}

\value{
list with following components
\item{u}{is a p x r matrix of kx row blocks uk (pk x r); uk'*uk = Identity}
\item{v}{is a q x r matrix of ky row blocks vh (qh x r); vh'*vh = Identity}
\item{s2}{is a kx x ky x r array; with r fixed, each matrix contains kxky values \eqn{(u_h'*x_{kh}*v_k)^2}, the partial (squared) singular values relative to \eqn{x_{kh}}}
}

\details{The first solution calculates kx+ky normed vectors: kx vectors
uk[:,1] of Rpk associated to ky vectors vh[,1]'s of Rqh, by maximizing
\eqn{\sum_k \sum_h (u_k[,1]'*x_{kh}*v_h[,1])^2}, with kx+ky norm
constraints.  A value \eqn{(u_k[,1]'*x_{kh}*v_h[,1])^2} measures the
relative link between \eqn{R^{p_k}} and \eqn{R^{q_h}} associated to the
block xkh.

The second solution is obtained from the same criterion, but after
replacing each xhk by xkh-xkh*vh*vh'-uk*uk'xkh+uk*uk'xkh*vh*vh'.  And
so on for the successive solutions 1,2,...,r .  The biggest number of
solutions may be r=inf(pk,qh), when the xkh's are supposed with full
rank; then rmax=min([min(K),min(H)]).

When K=p (or H=q, with t(x)), svdcp function is better.  When H=q and
K=p, it is the usual svd (with squared singular values).

Convergence of algorithm may be not global. So the below proposed
initialisation of the algorithm may be not very suitable for some data
sets.  Several different random initialisations with normed vectors
might be considered and the best result then choosen
}

\references{
Kissita G., Analyse canonique generalisee avec tableau de reference
generalisee. Thesis, Ceremade Paris 9 Dauphine (2003)
}

\examples{
x<-matrix(runif(200),10,20)
s2<-svdbip2(x,c(3,4,3),c(5,5,10),3);s2$s2
s1<-svdbip(x,c(3,4,3),c(5,5,10),3);s1$s2
}

