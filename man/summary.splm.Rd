\name{summary.splm}
\alias{summary.splm}
\title{summary method for class splm}

\description{
Method for summarizing the results of objects of class \code{'splm'}
}

\usage{

\method{summary}{splm}(object,...)
}

\arguments{
\item{object}{an object of class \code{'splm'}}
\item{...}{additional arguments to be passed}
}


\details{
The summary function \code{summary.splm} 
returns an objects of class 'splm' 
organized in a coefficient matrix.

Also a matrix for the error components, 
or the spatial coefficients will be generated 
depending on the estimated model. 

When the \code{'splm'} is produced by  
the function 'spsegm', the summary 
will be generated looping over the number
of equations in the system.
}

\value{
An object of class \code{'summary.splm'}

}

\seealso{
  \code{spml}, \code{spgm}
}
\author{Giovanni Millo, Gianfranco Piras}

\examples{
data(Produc, package = "plm") 
data(usaww)
GM <- spgm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc,
            listw=usaww, moments = "fullweights", spatial.error = TRUE)
summary(GM)
}

\keyword{spatial}
