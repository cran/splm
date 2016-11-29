\name{slmtest}
\alias{slmtest}
\alias{slmtest.formula}
\alias{slmtest.plm}

\title{Locally robust panel Lagrange Multiplier tests for spatial dependence}
\description{
  Locally robust LM tests for spatial lag (error) correlation sub
  spatial error (lag) correlation in panel models
}
\usage{
slmtest(x,...)
\method{slmtest}{formula}(formula, data, listw, model="pooling",
 test=c("lme","lml","rlme","rlml"), index=NULL, ...)
\method{slmtest}{plm}(x, listw,
 test=c("lme","lml","rlme","rlml"), ...)
}
\arguments{
\item{formula}{an object of class  \code{formula}}
\item{data}{a \code{data.frame} or \code{pdata.frame} containing the
    variables in the model}
\item{x}{an object of class \code{plm}}
\item{listw}{either a \code{matrix} or a \code{listw} representing the
  spatial structure}
\item{model}{a character value specifying the transformation to be
  applied to the data.} 
\item{test}{one of \code{c("lme","lml","rlme","rlml")}, the
  test to be performed.}
  \item{index}{either NULL (default) or a character vector to identify the indexes among the columns of the \code{data.frame}}
\item{...}{additional arguments to be passed}
}
\details{
  This tests are panel versions of the locally robust LM tests of
  Anselin et al. (1996), based on a pooling assumption: i.e., they do
  not allow for any kind of individual effect. Therefore it is advisable
  to employ a within transformation whenever individual effects cannot
  be ruled out.

  It must be kept in mind that these locally robust procedures have been
  designed for situations in which the "other" effect is not of
  substantial magnitude, and can behave suboptimally otherwise.

  Four tests are available to be chosen through the \code{test}
  argument: \code{"lml"} for "LM lag" and, respectively, \code{"lme"}
  for "LM error" are the standard, non-robust versions, obtained simply
  pooling the cross-sectional versions; \code{"rlml"} and \code{"rlme"}
  are, respectively, the locally robust test for lag, allowing for a
  spatial error; and for error, allowing for a spatial lag.
  
  The \code{model} argument, specified according to the standards of
\code{plm}, is passed on internally and employed to determine the panel
data transformation to be applied before calculating the test. Defaults
to \code{"pooling"} (no transformation).

  }
\value{
an object of class \code{htest}
}

\references{Anselin, L., Bera, A.K., Florax, R. and Yoon, M.J. (1996)
Simple diagnostic tests for spatial dependence.
\emph{Regional Science and Urban Economics}, \bold{26}(1), 77-104.
Elhorst, J.P. (2014) Spatial Panel data Models, in Spatial Econometrics
(Springer) 37-93.}

\author{Giovanni Millo}
%\seealso{\code{detest}}

\examples{
data(Produc, package="plm")
data(usaww)
fm <- log(gsp)~log(pcap)+log(pc)+log(emp)+unemp
## robust LM test for spatial error sub spatial lag
## model on original data, pooling hypothesis
slmtest(fm, data=Produc, listw = usaww, test="rlme")
## model on within-transformed (time-demeaned) data,
## eliminates individual effects
slmtest(fm, data=Produc, listw = usaww, test="rlme",
  model="within")
}

\keyword{htest}
