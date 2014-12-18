\name{spreml}
\alias{spreml}
\title{Spatial Panel Model with Random Effects by Maximum Likelihood}
\description{Maximum likelihood (ML) estimation of spatial panel models
  with random effects and serial error correlation.}


\usage{
spreml(formula, data, index = NULL, w, w2=w, lag = FALSE,
          errors = c("semsrre", "semsr", "srre", "semre",
          "re", "sr", "sem","ols", "sem2srre", "sem2re"),
          pvar = FALSE, hess = FALSE, quiet = TRUE,
          initval = c("zeros", "estimate"),
          x.tol = 1.5e-18, rel.tol = 1e-15, ...) 
}

\arguments{
  \item{formula}{a symbolic description of the model to be estimated}
  \item{data}{an object of class \code{data.frame} or \code{pdata.frame}. A data frame containing the variables
  in the model. When the object is a \code{data.frame}, the first two columns shall contain the indexes, unless otherwise specified. See \code{index}}
  \item{index}{if not NULL (default), a character vector to identify the indexes among the columns of the \code{data.frame}}
  \item{w}{an object of class \code{listw} or a \code{matrix}. It represents the spatial weights to be used in estimation.}
  \item{w2}{an object of class \code{listw} or a
  \code{matrix}. Second set of spatial weights for estimation, if
  different from the first (e.g., in a 'sarar' model).}
  \item{lag}{default=\code{FALSE}. If \code{TRUE}, a spatial lag of the dependent variable is added.}
  \item{errors}{Specifies the error covariance structure. See details.}
  \item{pvar}{legacy parameter here only for compatibility.}
  \item{hess}{default=\code{FALSE}. If \code{TRUE} estimate the
  covariance for beta_hat by numerical Hessian instead of GLS at optimal
  values.}
  \item{quiet}{default=\code{TRUE}. If \code{FALSE}, report function and
  parameters values during optimization.}
  \item{initval}{one of \code{c("zeros", "estimate")}, the initial values for
  the parameters. If \code{"zeros"} a vector of zeros is used. if
  \code{"estimate"} the initial values are retreived from the estimation
  of the nested specifications. Alternatively, a numeric vector can be
  specified.}
  \item{x.tol}{control parameter for tolerance. See \code{nlminb} for details.}
  \item{rel.tol}{control parameter for relative tolerance. See \code{nlminb} for details.}
  \item{...}{additional arguments to pass over to other functions, e.g. \code{method}.}
}
\details{Second-level wrapper for estimation of random effects models
  with serial and spatial correlation. The specifications without serial
  correlation (no "sr" in \code{errors}) can be called through
  \code{spml}, the extended ones only through \code{spreml}. 
  The models are estimated by two-step Maximum Likelihood.
  Abbreviations in \code{errors} correspond to: \code{"sem"}
  Anselin-Baltagi type spatial autoregressive error: if
  present, random effects are not spatially correlated;  \code{"sem2"}
  Kapoor, Kelejian and Prucha-type spatial autoregressive error model
  with spatially correlated random effects; \code{"sr"} serially
  correlated remainder errors; \code{"re"} random effects; \code{"ols"}
  spherical errors (usually combined with \code{lag=T}).
  The optimization \code{method} can be passed on as optional
  parameter. Default is \code{"nlminb"}; all constrained optimization
  methods from \code{maxLik} are allowed (\code{"BFGS", "NM", "SANN"})
  but the latter two are still experimental.
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
  Millo, G. (2014)
  Maximum likelihood estimation of spatially and serially correlated
  panels with random effects. \emph{Computational Statistics and Data
  Analysis}, \bold{71}, 914--933.
}

\author{Giovanni Millo}

\seealso{\code{spml}}

\examples{
data(Produc, package = "plm")
data(usaww)
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
## random effects panel with spatial lag and serial error correlation
## optimization method set to "BFGS"
sarsrmod <- spreml(fm, data = Produc, w = usaww, errors="sr", lag=TRUE,
                    method="BFGS")
summary(sarsrmod)
}

\keyword{spatial}
