\name{effects.splm}
\alias{effects.splm}
\alias{print.effects.splm}
\title{method for extracting fixed effects}

\description{
Methods used for extracting fixed effects from objects of class \code{splm} where \code{type} is one of 
"fixed effects lag" or "fixed effects error"
}

\usage{
\method{effects}{splm}(object,...)
}

\arguments{
\item{object}{an object of class \code{'splm'}}
\item{...}{additional arguments to be passed over}
}

\details{
If the argument \code{object} is not of class \code{splm} 
the function will terminate with an error.

If the argument \code{object} is of class \code{splm} 
but \code{type} is not one of 
"fixed effects lag" or "fixed effects error",
the function will terminate with an error.

}


\value{
  An object of class \code{effects.splm}
  \item{res}{a list whose elements are various type of fixed effects and the intercept (when present)}
}

\references{
  Elhorst, J.P. (2003) 
Specification and estimation of spatial panel data models,
    \emph{International Regional Science Review}, \bold{26}, pages 244--268.

  Elhorst, J.P. (2009) 
Spatial panel data models,
    \emph{In} Fischer, M.M. and Getis, A. (eds),
    \emph{Handbook of Applied Spatial Analysis} Springer, Berlin.

}

\author{Gianfranco Piras}

\seealso{
\code{spml} \code{summary.effects.splm}
}
\examples{
data(Produc, package = "plm")
data(usaww)
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
err <- spml(fm, data = Produc, listw = spdep::mat2listw(usaww), model="within")
summary(err)
eff <- effects(err) 
print(eff)
}
\keyword{spatial}
