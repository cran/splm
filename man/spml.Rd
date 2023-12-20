\name{spml}
\alias{spml}
\alias{impacts}
\alias{impacts.splm_ML}
\alias{impacts.splm_GM}
\title{Spatial Panel Model by Maximum Likelihood}
\description{Maximum likelihood (ML) estimation of spatial panel models, possibly with fixed or random effects.}


\usage{
spml(formula, data, index=NULL, listw, listw2=listw, na.action,
                 model=c("within","random","pooling"),
                 effect=c("individual","time","twoways"),
                 lag=FALSE, spatial.error=c("b","kkp","none"),
                 ...)
\method{impacts}{splm_ML}(obj, listw = NULL,
                         time = NULL, ...,
                         tr = NULL, R = 200,
                         type = "mult",
                         empirical = FALSE, Q = NULL)
\method{impacts}{splm_GM}(obj, ..., tr=NULL, 
                            R=NULL, listw=NULL,
                            type = "mult",
                            time = NULL,
                            evalues=NULL, tol=1e-6, 
                            empirical=FALSE, Q=NULL, 
                            KPformula = FALSE, prt = TRUE)
}

\arguments{
  \item{formula}{a symbolic description of the model to be estimated}
  \item{data}{an object of class \code{data.frame} or \code{pdata.frame}. A data frame containing the variables
  in the model. When the object is a \code{data.frame}, the first two columns shall contain the indexes, unless otherwise specified. See \code{index}}
  \item{index}{if not NULL (default), a character vector to identify the indexes among the columns of the \code{data.frame}}
  \item{listw}{an object of class \code{listw} or a \code{matrix}. It represents the spatial weights to be used in estimation.}
  \item{listw2}{an object of class \code{listw} or a
  \code{matrix}. Second of set spatial weights for estimation, if
  different from the first (e.g., in a 'sarar' model).}
  \item{na.action}{see \pkg{spdep} for more details.}
  \item{model}{one of \code{c("within", "random", "pooling").}}
  \item{effect}{one of \code{c("individual","time","twoways")}; the
  effects introduced in the model.}
  \item{lag}{default=\code{FALSE}. If \code{TRUE}, a spatial lag of the dependent variable is added.}
  \item{spatial.error}{one of \code{c("b","kkp","none")}. The type of
  spatial error in the specification, if any. See details.}
  \item{...}{additional argument to pass over to other functions}
  \item{obj}{fitted model object}
  \item{time}{??time??}
  \item{tr}{A vector of traces of powers of the spatial weights matrix created using 'trW', for approximate impact measures}
  \item{R}{If given, simulations are used to compute distributions for the impact measures, returned as 'mcmc' objects}
  \item{type}{Either "mult" (default) for powering a sparse matrix (with moderate or larger N, the matrix becomes dense, and may lead to swapping), or "MC" for Monte Carlo simulation of the traces (the first two simulated traces are replaced by their analytical equivalents), or "moments" to use the looping space saving algorithm proposed by Smirnov and Anselin (2009) - for "moments", 'W' must be symmetric, for row-standardised weights through a similarity transformation}
  \item{empirical}{Argument passed to 'mvrnorm' (default FALSE)}
  \item{Q}{default NULL, else an integer number of cumulative power series impacts to calculate if 'tr' is given}
  \item{evalues}{vector of eigenvalues of spatial weights matrix for impacts calculations}
  \item{tol}{Argument passed to 'mvrnorm'}
  \item{KPformula}{not yet implemented}
  \item{prt}{not yet implemented}
}
\details{
  The models are estimated by two-step Maximum Likelihood.
  Further optional parameters to be passed on to the estimator may be:
pvar: if \code{TRUE} the \code{pvar} function is called
hess: if \code{TRUE} use numerical Hessian instead of GLS for the
standard errors of the estimates
quiet: if \code{FALSE} report function and parameters values during
optimization
initval: one of \code{c("zeros", "estimate")}, the initial values for
  the parameters. If \code{"zeros"} a vector of zeros is used. if
  \code{"estimate"} the initial values are retreived from the estimation
  of the nested specifications. Alternatively, a numeric vector can be
  specified.
  x.tol: Tolerance. See \code{nlminb} for details.
  rel.tol: Relative tolerance. See \code{nlminb} for details.
}
\value{
  An object of class \code{"splm"}.
  \item{coefficients}{coefficients estimate of the model parameters }
  \item{arcoef}{the coefficient for the spatial lag on \code{y}}
  \item{errcomp}{the estimates of the error variance components}
  \item{vcov}{the asymptotic variance covariance matrix of the estimated coefficients}
  \item{vcov.arcoef}{the asymptotic variance of the
  estimated spatial lag parameter}
  \item{vcov.errcomp}{the asymptotic variance covariance matrix of the
  estimated error covariance parameters}
  \item{type}{'random effects ML'}
  \item{residuals}{the model residuals}
  \item{fitted.values}{the fitted values, calculated as \eqn{\hat{y}=X \hat{\beta}}}
  \item{sigma2}{GLS residuals variance}
  \item{model}{the matrix of the data used}
  \item{call}{the call used to create the object}
  \item{logLik}{the value of the log likelihood function at the optimum }
  \item{errors}{the value of the \code{errors} argument}
}
\references{
Baltagi, B.H., Song, S.H., Jung B. and Koh, W. (2007)
Testing panel data regression models with spatial and serial error correlation.
\emph{Journal of Econometrics}, \bold{140}, 5-51.

Millo, G., Piras, G. (2012)
splm: Spatial Panel Data Models in R.
\emph{Journal of Statistical Software}, \bold{47(1)}, 1--38.
URL http://www.jstatsoft.org/v47/i01/.
}

\author{Giovanni Millo}

\seealso{\code{spgm}}

\examples{
data(Produc, package = "plm")
data(usaww)
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
## the two standard specifications (SEM and SAR) one with FE
## and the other with RE:
## fixed effects panel with spatial errors
fespaterr <- spml(fm, data = Produc, listw = spdep::mat2listw(usaww),
                   model="within", spatial.error="b", Hess = FALSE)
summary(fespaterr)
## random effects panel with spatial lag
respatlag <- spml(fm, data = Produc, listw = spdep::mat2listw(usaww),
                   model="random", spatial.error="none", lag=TRUE)
summary(respatlag)
## calculate impact measures
#impac1 <- impacts.splm(respatlag, listw = spdep::mat2listw(usaww, #style = "W"), time = 17)
#summary(impac1, zstats=TRUE, short=TRUE)
}

\keyword{spatial}
