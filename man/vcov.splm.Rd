\name{vcov.splm}
\alias{vcov.splm}

\title{Covariance extractor method for \code{splm} objects}
\description{
  Covariance extractor method for \code{splm} objects. Seldom used as
  such but needed, e.g., for
  interoperability with testing functions in lmtest and car.
}
\usage{
\method{vcov}{splm}(object, ...)
}
\arguments{
\item{object}{an object of class  \code{splm}}
\item{...}{additional arguments to be passed; currently not used}
}

\value{
a covariance matrix of beta coefficients
}


\author{Giovanni Millo}

\references{Zeileis, A. (2006) Object-Oriented Computation of Sandwich
  Estimators. \emph{Journal of Statistical Software}, \bold{16(9)}, 1-16.}

\examples{
## not run:
## data(Produc, package="plm")
## data(usaww)
## fm <- log(gsp)~log(pcap)+log(pc)+log(emp)+unemp
## sarremod <- spml(fm, data=Produc, listw = mat2listw(usaww),
##  model="random", lag=TRUE, spatial.error="none")
## ## compact representation of betas
## library(lmtest)
## coeftest(sarremod)
## ## linear hypothesis test
## library(car)
## lht(sarremod, "log(pcap)=log(pc)")
}

\keyword{vcov}
