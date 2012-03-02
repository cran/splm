\name{write.effects.splm}
\alias{write.effects.splm}
\title{method for writing a table with fixed effects}

\description{
Methods used for writing fixed effects to \code{filename}
}

\usage{
\method{write}{effects.splm}(x, file = "effects", ncolumns = if (is.character(x)) 1 else 5, 
    append = FALSE, sep = ",")
}

\arguments{
\item{x}{an object of class \code{'effects.splm'}}
\item{file}{name of the file, default set to \code{"effects"}}
\item{ncolumns}{the number of columns to write the data in}
\item{append}{if \code{TRUE} the data are appended to the connection}
\item{sep}{string to separate columns; default is = ","}
}

\details{
If the argument \code{object} is not of class \code{effects.splm}
the function will terminate with an error message.

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

\author{ Gianfranco Piras \email{gpiras@mac.com}}

\seealso{
\code{\link{spml}}
}
\examples{
data(Produc, package = "Ecdat")
data(usaww)
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
lag <- spml(fm, data = Produc, listw = mat2listw(usaww), model="within", effect="twoways", lag=TRUE, method="eigen")
summary(lag)
eff <- effects(lag) 
print(eff)
err <- spml(fm, data = Produc, listw = mat2listw(usaww), model="within", effect="time", spatial.error="b", method="eigen")
summary(err)
eff <- effects(err) 
write.effects.splm(eff)
}

\keyword{spatial}
