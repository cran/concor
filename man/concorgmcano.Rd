\name{concorgmcano}
\title{Canonical analysis of subsets Yj with subsets Xi}
\description{Canonical analysis of subsets Yj with subsets Xi. Relative valuations by squared correlations of the proximities of subsets Xi with subsets Yj. SUCCESSIVE SOLUTIONS}
\keyword{}
\keyword{}
\alias{concorgmcano}

\usage{
concorgmcano(x,px,y,py,r)
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
\item{rho2}{is a kx x ky x r array; for a fixed solution k, rho2[,,k]
  contains kxky squared correlations \eqn{rho2(cx[n*(i-1)+1:n*i,k],cy[n*(j-1)+1:n*j,k])}, simultaneously calculated between all the yj with all the xi}
}

\details{
For the first solution, \eqn{sum_i sum_j \mbox{rho2}(cx_i[,1],cy_j[,1])} is the
optimized criterion. The other solutions are calculated from the same
criterion, but with orthogonalities for having two by two zero
correlated the canonical components defined for each xi, and also for
those defined for each yj.  Each solution associates kx canonical
components to ky canonical components.  When kx =1 (px=p), take
concorcano function

This function uses the concorgm function
}

\references{
Kissita G., Analyse canonique generalisee avec tableau de reference generalisee. Thesis, Ceremade Paris 9 Dauphine (2003).
}

\examples{
x<-matrix(runif(50),10,5);y<-matrix(runif(90),10,9)
x<-scale(x);y<-scale(y)
cc<-concorgmcano(x,c(2,3),y,c(3,2,4),2)
diag(t(cc$cx[1:10,])\%*\%cc$cy[1:10,]/10)^2
cc$rho2[1,1,] 
}

