\name{concoreg}
\title{Redundancy of sets yj by one set x}
\description{Regression of several subsets of variables Yj by another set X. SUCCESSIVE SOLUTIONS}
\keyword{}
\keyword{}
\alias{concoreg}

\usage{
concoreg(x,y,py,r)
}

\arguments{
\item{x}{is a \eqn{n \times p}{n x p} matrix of p centered explanatory variables}
\item{y}{is a \eqn{n \times q}{n x q} matrix of q centered variables}
\item{py}{is a row vector which contains the numbers \eqn{q_i,
    i=1,...,ky}, of the ky subsets \eqn{y_i} of y : \eqn{\sum_i q_i} =
    sum(py) = q. py is the partition vector of y} \item{r}{is the wanted
    number of successive solutions} }

\value{
list with following components
\item{cx}{the \eqn{n \times r}{n x r} matrix of the r explanatory components}
\item{v}{is a \eqn{q \times r}{q x r} matrix of ky row blocks \eqn{v_i} (\eqn{q_i \times r}{q_i x r}) of axes in Rqi relative to yi; \eqn{v_i'*v_i = \mbox{Id}}}
\item{V}{is a \eqn{q \times r}{q x r} matrix of axes in Rq relative to y; \eqn{V'*V = \mbox{Id}}}
\item{varexp}{is a \eqn{ky \times r}{ky x r} matrix; each column k contains ky explained variances \eqn{\rho(cx[,k],y_i*v_i[,k])^2 \mbox{var}(y_i*v_i[,k])}}
}

\details{

The first solution calculates 1+ky normed vectors: the component
cx[,1] in \eqn{R^n} associated to the ky vectors vi[,1]'s of \eqn{R^{q_i}}, by
maximizing \eqn{varexp1=\sum_i \rho(cx[,1],y_i*v_i[,1])^2 \mbox{var}(y_i*v_i[,1]))}, with
\eqn{1+ky} norm constraints. A explanatory component cx[,k] is associated to
ky partial explained components yi*vi[,k] and also to a global explained
component y*V[,k]. \eqn{\rho(cx[,k],y*V[,k])^2 \mbox{var}(y*V[,k])= \mbox{varexpk}}.  The
total explained variance by the first solution is maximal.

The second solution is obtained from the same criterion, but after
replacing each yi by \eqn{y_i-y_i*v_i[,1]*v_i[,1]'}.  And so on for the
successive solutions 1,2,...,r .  The biggest number of solutions may
be \eqn{r=inf(n,p,q_i)}, when the matrices x'*yi are supposed with full rank.
For a set of r solutions, the matrix (cx)'*y*V is diagonal : "on
average", the explanatory component of one solution is only linked
with the components explained by this explanatory, and is not linked
with the explained components of the other solutions.  The matrices
\eqn{(cx)'*y_j*v_j} are triangular : the explanatory component of one solution
is not linked with each of the partial components explained in the
following solutions.  The definition of the explanatory components
depends on the partition vector py from the second solution.

This function is using concor function 
}

\references{
Hanafi & Lafosse (2001) Generalisation de la regression lineaire
simple pour analyser la dependance de K ensembles de variables avec un
K+1 eme. Revue de Statistique Appliquee vol.49, n.1.

Chessel D. & Hanafi M. (1996) Analyses de la Co-inertie de K nuages de
points.  Revue de Statistique Appliquee vol.44, n.2. (this ACOM
analysis of one multiset is obtained by the command : concoreg(Y,Y,py,r))
}

\examples{
x<-matrix(runif(50),10,5);y<-matrix(runif(90),10,9)
x<-scale(x);y<-scale(y)
co<-concoreg(x,y,c(3,2,4),2)
((t(co$cx[,1])\%*\%y[,1:3]\%*\%co$v[1:3,1])/10)^2;co$varexp[1,1]
t(co$cx)\%*\%co$cx /10
diag(t(co$cx)\%*\%y\%*\%co$V/10)^2
sum(co$varexp[,1]);sum(co$varexp[,2])
}

