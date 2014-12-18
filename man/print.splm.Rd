\name{print.splm}
\alias{print.splm}
\title{print method for class splm}

\description{
Method to print objects of class \code{summary.splm} and \code{splm}
}

\usage{
\method{print}{splm}(x, digits = max(3,getOption("digits") -3), ...)
}

\arguments{
\item{x}{an object of class \code{splm}}
\item{digits}{minimal number of significant digits, see \code{print.default}}
\item{...}{additional arguments to be passed}
}


\details{
The summary function \code{summary.splm} 
returns an objects of class 'splm' 
organized in a coefficient matrix.

Also a matrix for the error components, 
or the spatial coefficients will be generated 
depending on the estimated model. 
}

\seealso{
  \code{spml}, \code{spgm}
}
\author{Giovanni Millo, Gianfranco Piras}

\examples{
data(Produc, package = "plm") 
data(usaww)
spremod<-spml(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc,
listw = mat2listw(usaww), model="random", lag=TRUE, spatial.error="none")
summary(spremod)
}
\keyword{spatial}
