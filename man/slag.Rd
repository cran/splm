\name{slag}
\alias{slag}
\alias{slag.default}
\alias{slag.pseries}

\title{Spatial lag operator}
\description{
  Spatial lagging method for vectors or \code{pseries} objects.
}
\usage{
\method{slag}{pseries}(x, listw, maxlag, ...)
}
\arguments{
  \item{x}{an object of class  \code{pseries}}
  \item{listw}{an object of class \code{listw}}
  \item{maxlag}{the spatial lag order (including lower)} 
\item{...}{additional arguments to be passed}
}

\value{
a \code{pseries}
}


\author{Giovanni Millo}

\examples{
data(Produc, package="plm")
data(usaww)
usalw <- spdep::mat2listw(usaww)
fm <- log(gsp)~log(pcap)+log(pc)+log(emp)+unemp+slag(log(pcap),
  listw=usalw)
slxmod <- spreml(fm, data=Produc, w = usaww,
  model="pooling", lag=FALSE, errors="ols")
}

\keyword{slag}
