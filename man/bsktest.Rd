\name{bsktest}
\alias{bsktest}
\alias{bsktest.formula}

\title{Baltagi, Song and Koh LM test for spatial panels}
\description{
  Baltagi, Song and Koh marginal or conditional LM test for
   spatial error correlation or random effects in panel models
}
\usage{
bsktest(x,...)
\method{bsktest}{formula}(x, data, index=NULL, listw,
test=c("LMH","LM1","LM2","CLMlambda","CLMmu"), 
standardize=TRUE, ...)
}
\arguments{
\item{x}{a \code{formula}}
  \item{data}{a \code{data.frame} or \code{pdata.frame} containing the variables
  in the model}
  \item{index}{either NULL (default) or a character vector to identify the indexes among the columns of the \code{data.frame}}
\item{listw}{a \code{listw} representing the spatial structure}
\item{test}{one of \code{c("LMH","LM1","LM2","CLMlambda","CLMmu")}, the
  test to be performed}
\item{standardize}{whether to standardize the test statistic or not (applies only to LM1 and LM2)}
\item{...}{additional arguments to be passed}
}

\value{
an object of class \code{htest}
}
\references{Baltagi, B.H., Song, S.H. and Koh, W. (2003)
Testing panel data regression models with spatial error correlation.
\emph{Journal of Econometrics}, \bold{117}, 123--150.

Millo, G., Piras, G. (2012)
splm: Spatial Panel Data Models in R.
\emph{Journal of Statistical Software}, \bold{47(1)}, 1--38.
URL http://www.jstatsoft.org/v47/i01/.}
\author{Gianfranco Piras}
\seealso{\code{sphtest}}

\examples{
data(Produc, package="Ecdat")
data(usaww)
fm <- log(gsp)~log(pcap)+log(pc)+log(emp)+unemp
bsktest(fm,data=Produc, listw = mat2listw(usaww),
  test="LM1", standardize=TRUE)
}

\keyword{htest}
