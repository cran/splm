\name{bsjktest}
\alias{bsjktest}
\alias{bsjktest.formula}

\title{Baltagi, Song, Jung and Koh LM test for spatial panels}
\description{
  Baltagi, Song, Jung and Koh joint or conditional LM test for
   spatial error correlation or serial correlation sub spatial, serial correlation and random effects in panel models
}
\usage{
bsjktest(x,...)
\method{bsjktest}{formula}(x, data, index=NULL, listw,
test=c("C.1","C.2","C.3","J"), ...)
}
\arguments{
\item{x}{an object of class  \code{formula}}
  \item{data}{a \code{data.frame} or \code{pdata.frame} containing the variables in the model}
  \item{index}{either NULL (default) or a character vector to identify the indexes among the columns of the \code{data.frame}}
\item{listw}{either a \code{matrix} or a \code{listw} representing the spatial structure}
\item{test}{one of \code{c("C.1","C.2","C.3","J")}, the
  test to be performed. "C.3" is not implemented yet.}
\item{...}{additional arguments to be passed}
}

\value{
an object of class \code{htest}
}
\references{Baltagi, B.H., Song, S.H., Jung B. and Koh, W. (2007)
Testing panel data regression models with spatial and serial error correlation.
\emph{Journal of Econometrics}, \bold{140}, 5-51.}
\author{Giovanni Millo}
\seealso{\code{bsktest}}

\examples{
data(Produc, package="plm")
data(usaww)
fm <- log(gsp)~log(pcap)+log(pc)+log(emp)+unemp
bsjktest(fm, data=Produc, listw = usaww, test="C.1")
}

\keyword{htest}
