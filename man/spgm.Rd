\name{spgm}
\alias{spgm}
\alias{ivsplm} 
\alias{spsarargm} 
\alias{sperrorgm}
\title{GM estimation of spatial panel data models}

\description{GM estimation of panel data models with spatially correlated 
errors components of the form: 

\deqn{ y_N(t) = \lambda W y + X_N(t) \beta + u_N(t) }

\deqn{ u_N(t) = \rho W_N u_N(t) + \epsilon(t)}

\deqn{ \epsilon_N = (e_T \otimes I_N ) \mu_N + \nu_N }

where \eqn{ \rho}, and the variance components \eqn{\sigma^2_\mu} and \eqn{\sigma^2_\nu} 
are estimated by GM, and the model coefficients by a feasible GLS estimator. The model can also include 
additional (other than the spatial lag) endogenous variables. 
}

\usage{
spgm(formula, data=list(), index=NULL, listw =NULL, listw2 = NULL,
         model=c("within","random"), lag = FALSE, spatial.error=TRUE,
         moments = c("initial", "weights", "fullweights"), endog = NULL, 
         instruments= NULL, lag.instruments = FALSE, verbose = FALSE, 
         method = c("w2sls", "b2sls", "g2sls", "ec2sls"), control = list(), 
         optim.method = "nlminb", pars = NULL)
}

\arguments{
  \item{formula}{a description of the model to be fit. The details of model specification are given
  for \code{lm} }
  \item{data}{an object of class \code{data.frame} or \code{pdata.frame}. An optional data frame containing the variables
  in the model. When the obect is a \code{data.frame}, the first two columns may contain the indexes. See \code{index}}
  \item{index}{if not NULL (default), a character vector to identify the indexes among the columns of the \code{data.frame}}
  \item{listw}{an object of class \code{listw}, \code{matrix}, or \code{Matrix}}. 
  \item{listw2}{an object of class \code{listw}, \code{matrix}, or \code{Matrix}. Only 
  if both \code{lag} and \code{spatial.error} are both \code{TRUE}}
  \item{model}{One of \code{"within"} or \code{"random"}. The assumption
    made on the individual effects}
  \item{lag}{if \code{TRUE} a spatial lag of the dependent variable is
    added to the regression equation}
  \item{spatial.error}{a logic vector. If \code{TRUE} the spatial autoregressive error term is added to the model and an estimate for \eqn{\rho} is produced}
  \item{moments}{\code{"initial"} (default) defines the set of GM estimator to be used. Alternatives are  
  \code{"weights"} and \code{"fullweights"} (See Details)}
  \item{endog}{additional endogenous variables. Default \code{NULL}. If not \code{NULL} should be specified 
as a formula with no dependent variable (endog = ~ x1 + x2). Note the ~ before the expression. }
  \item{instruments}{external instruments. Default \code{NULL}. If not \code{NULL} should be specified 
as a formula with no dependent variable (instruments = ~ x1 + x2). Note the ~ before the expression.}
\item{lag.instruments}{should the external instruments be spatially lagged?}
  \item{verbose}{default \code{FALSE}, If \code{TRUE} reports function values during optimization}  
  \item{method}{One of \code{"w2sls"}, \code{"b2sls"}, \code{"g2sls"}, \code{"ec2sls"}. (See Details)}
  \item{control}{a list of control parameters for the optimization}
    \item{optim.method}{default set to \code{"nlminb"}. or optionally a method passed to \code{optim} to use an alternative optimizer.}
    \item{pars}{initial values of the parameter \code{rho} and \code{sigmav}. The default for \code{rho} is to 
    start from a regression of the spatially lagged residuals on the residuals (depending on the model). 
    For \code{sigmav} the starting value is the variance of the residuals (again this depends on the model).  }
}


\details{ 
The function is a very general interface to estimate various nested specifications of the 
general model including additional endogenous variables described above. 
When both \code{spatial.error} and \code{lag} are \code{FALSE} the model reduces to a 
panel data model with an additional endogeneous variable. The function then uses \code{ivsplm} 
to perform the Instrumental Variables and two-stage least squares for panel data model.  \code{method = "w2sls"} 
corresponds to the fixed effects estimator, \code{method = "b2sls"} to the between effects model,
\code{method = "g2sls"} to the GLS random effects model, and \code{method = "ec2sls"} to teh Baltagi's 
EC2SLS. 

When \code{spatial.error} is \code{TRUE} and \code{lag} is \code{FALSE} the model is one with spatially autocorrelated
error components. If \code{effects} is \code{"random"}, the Kapoor et al. (2007) GM estimator is performed and
the residuals in the first step come from an OLS regression. 
When \code{moments} is \code{"initial"}, the initial estimator is calculated. This first set 
of GM estimators is based only on a subset of the moments conditions and assigns 
equal weigths to each of them. 
When \code{moments} is \code{"fullweights"}, the second set of GM estimators is calculated. This estimator 
is based on the full set of moments conditions. It also involves the expression for the variance 
covariance matrix of the sample moments calculated under the assumption of 
normally distributed innovations. The calculation of the trace terms in the expression 
of the variance covariance matrix of the sample moments 
uses codes from the \code{Matrix} package. 
When \code{moments} is\code{"weights"}, the third set of GM estimator is used. This is motivated by computational 
issues. The procedure is analogous to the second one but uses a simplified expression for the 
variance covariance matrix of the sample moments.
If \code{effects} is \code{"fixed"}, the initial 
estimator is a within estimator and the moments conditions of Kapoor et al. (2007) are modified accordingly.

Finally, when both \code{spatial.error} and \code{lag} are \code{TRUE} the complete model is estimated (with or without 
additional endogenous variables). OLS residuals are no longer consistent because of the spatially lagged dependent variable.
If \code{effects} is \code{"random"}, two initial estimators are computed: a within two-stage least squares
and a between two stage least squares. The two sets of corresponding residuals are used in the spatial generalized moments estimator
(GM) where the moments conditions of Kapoor et al. (2007) are again modified accordingly.
If \code{effects} is \code{"fixed"}, the initial 
estimator is a within two stage least squares estimator and the moments conditions of Kapoor et al. (2007) are modified accordingly.

Note that for the random effects models, \eqn{\sigma^2_\mu} is not reported. \eqn{\sigma^2_1} is reported instead. 
However, a value for \eqn{\sigma^2_\mu} can easily be obtained from:
\deqn{\sigma^2_1 = \sigma^2_\nu + T \sigma^2_\mu}
The function also produces an estimate for \eqn{\theta} which is a
function of the variance components.
}


\value{
  An object of class \code{"splm"}.
  \item{coefficients}{GLS coefficients estimate of the model parameters}
  \item{vcov}{the variance covariance matrix of the estimated coefficients}
  \item{residuals}{the GLS residuals}
  \item{fitted.values}{difference between response variable and residuals}
  \item{sigma2}{GLS residuals variance}
  \item{type}{'random effect GM'}
  \item{rho}{a vector including the spatial parameter and the variance components (see Details)}  
  \item{model}{the matrix of the data used}
  \item{call}{the call used to create the object}
}

\references{
  Kapoor, M., Kelejian, H.H. and Prucha, I.R. (2007) 
  Panel data model with spatially correlated error components,
    \emph{Journal of Econometrics}, \bold{140}, pages 97--130. 

  Mutl, J., and Pfaffermayr, M. (2011) 
  The Hausman test in a Cliff and Ord panel model,
    \emph{Econometrics Journal}, \bold{14}, pages 48--76. 
    
  Kelejian, H.H. and Prucha, I.R. (1999) 
A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model,
    \emph{International Economic Review}, \bold{40}, pages 509--533.
    
      Kelejian, H.H. and Prucha, I.R. (1999) 
A Generalized Spatial Two Stage Least Square Procedure for Estimating a Spatial Autoregressive
Model with Autoregressive Disturbances,
    \emph{Journal of Real Estate Finance and Economics}, \bold{17},
    pages 99--121.

    Millo, G., Piras, G. (2012)
splm: Spatial Panel Data Models in R.
\emph{Journal of Statistical Software}, \bold{47(1)}, 1--38.
URL http://www.jstatsoft.org/v47/i01/.
}

\author{Gianfranco Piras}


\examples{
data(Produc, package = "plm") 
data(usaww)
GM <- spgm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc,
           listw = usaww, moments="fullweights", spatial.error = TRUE)
summary(GM)
}

\keyword{spatial}
