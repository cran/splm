\name{splm-package}
\alias{splm-package}
\alias{splm}
\docType{package}
\title{
Spatial panel models: estimation and testing
}
\description{
A comprehensive toolset for ML and GM estimation and diagnostic testing of econometric models for spatial panel data.
}
\details{
\tabular{ll}{
Package: \tab splm\cr
Type: \tab Package\cr
Version: \tab 0.2-1\cr
Date: \tab 2010-06-08\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}

\author{
Giovanni Millo and Gianfranco Piras

Maintainer: Giovanni Millo <giovanni.millo@generali.com>
}

\references{
Anselin, L. and Le Gallo, J. and Jayet, H. (2008)
Spatial Panel Econometrics, \emph{In}
Matyas, L. and Sevestre, P. (eds), \emph{The econometrics of Panel
Data, Fundamentals and Recent Developments in Theory and Practice (3rd
Edition)}, pages 624--660. Springer-Verlag, Berlin Heidelberg.

Baltagi, B.H., Song, S.H., Jung B. and Koh, W. (2007)
Testing panel data regression models with spatial and serial error correlation.
\emph{Journal of Econometrics}, \bold{140}, 5-51.

Baltagi, B.H., Song, S.H. and Koh, W. (2003)
Testing panel data regression models with spatial error correlation.
\emph{Journal of Econometrics}, \bold{117}, 123--150.

  Elhorst, J.P. (2003)
Specification and estimation of spatial panel data models,
    \emph{International Regional Science Review}, \bold{26}, pages 244--268.

  Elhorst, J.P. (2009)
Spatial panel data models,
    \emph{In} Fischer, M.M. and Getis, A. (eds),
    \emph{Handbook of Applied Spatial Analysis} Springer, Berlin.

  Kapoor, M., Kelejian, H.H. and Prucha, I.R. (2007) 
  Panel data model with spatially correlated error components,
    \emph{Journal of Econometrics}, \bold{140}, pages 97--130. 
    
  Kelejian, H.H. and Prucha, I.R. (1999) 
A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model,
    \emph{International Economic Review}, \bold{40}, pages 509--533.
    
      Kelejian, H.H. and Prucha, I.R. (1999) 
A Generalized Spatial Two Stage Least Square Procedure for Estimating a Spatial Autoregressive
Model with Autoregressive Disturbances,
    \emph{Journal of Real Estate Finance and Economics}, \bold{17}, pages 99--121.

  Kelejian, H.H. and Prucha, I.R. (2004) 
Estimation of Simultaneous systems of spatially interrelated cross sectional equations,
\emph{Journal of Econometrics}, \bold{118}, pages 27--50.

Millo, G., Piras, G. (2012)
splm: Spatial Panel Data Models in R.
\emph{Journal of Statistical Software}, \bold{47(1)}, 1--38.
URL http://www.jstatsoft.org/v47/i01/.
}


\keyword{ package }

\examples{
data(Produc, package = "Ecdat") 
data(usaww)
Produc <- Produc[Produc$year<1975, ] 
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
GM<-spgm(log(gsp)~log(pcap)+log(pc)+log(emp) + unemp, data = Produc, listw = usaww, moments = "fullweights", spatial.error = TRUE)
summary(GM)
respaterr <- spml(fm, data = Produc, listw = mat2listw(usaww), model="random", spatial.error="b")
summary(respaterr)

}
