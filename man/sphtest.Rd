\name{sphtest}
\alias{sphtest}
\alias{sphtest.formula}
\alias{sphtest.splm}

\title{Hausman test for spatial panel data models}
\description{
Hausman specification test for spatial panel data models
}
\usage{
sphtest(x, ...)
\method{sphtest}{formula}(x, data, index = NULL, listw, spatial.model = c("lag", "error", "sarar"), method = c("ML", "GM"), errors = c("KKP", "BSK"),...)
\method{sphtest}{splm}(x, x2, ...)
}

\arguments{
\item{x}{an object of class  \code{formula} or \code{splm} }
\item{x2}{an object of class  \code{splm} }
  \item{data}{an object of class \code{\link{data.frame}} or \code{pdata.frame}. An optional data frame containing the variables
  in the model. When the obect is a \code{\link{data.frame}}, the first two columns may contain the indexes. See \code{index}}
\item{index}{if not NULL (default), a character vector to identify the indexes among the columns of the \code{\link{data.frame}}}
\item{listw}{an object of class \code{listw} created for example by \code{nb2listw} }
\item{spatial.model}{one of \code{c("lag","error","sarar")}, the model to be estimated (only lag, only error, both lag and error dependence)}
\item{method}{one of \code{c("ML","GM"}}
\item{errors}{one of \code{c("BSK","KPP"}. When \code{method} is \code{"ML"} defines the specification of the innovations}
\item{...}{additional arguments to be passed}
}

\value{
an object of class \code{htest}
}

\references{
  Millo, G., Piras, G. (2012)
splm: Spatial Panel Data Models in R.
\emph{Journal of Statistical Software}, \bold{47(1)}, 1--38.
URL http://www.jstatsoft.org/v47/i01/.}

\author{Gianfranco Piras}

\seealso{\code{\link{spgm}}}

\examples{
data(Produc, package="Ecdat")
data(usaww)
fm <- log(gsp)~log(pcap)+log(pc)+log(emp)+unemp
test1<-sphtest(fm,data=Produc, listw = mat2listw(usaww), spatial.model = "error", method="GM")
test1
mod1<- spgm(fm, data=Produc, listw = usaww, model = "random", spatial.error = TRUE, moments="fullweights")
mod2<- spgm(fm, data=Produc, listw = usaww, model = "within", spatial.error = TRUE)
test2<-sphtest(mod1, mod2)
test2
}
\keyword{htest}