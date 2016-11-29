\name{rwtest}
\alias{rwtest}
\alias{rwtest.formula}
\alias{rwtest.panelmodel}
\alias{rwtest.pseries}
\title{Randomization-based test of spatial dependence for panel models}

\description{
 Randomization-based test of spatial dependence for panel models, robust
 to global dependence induced by common factors and to persistence
 (serial correlation) in the data
}

\usage{
rwtest(x, \dots)
\method{rwtest}{formula}(x, data, w, index = NULL, model = NULL,
                         replications = 99, seed=NULL, order=1,
                         mc=1, test = c("rho", "cd", "sclm"),
                         alternative=c("twosided", "onesided",
                                       "symmetric"), \dots)
\method{rwtest}{panelmodel}(x, w, replications = 99, seed=NULL,
                         order=1, mc=1,
                         test = c("rho", "cd", "sclm"),
                         alternative=c("twosided", "onesided",
                                       "symmetric"), \dots)
\method{rwtest}{pseries}(x, w, replications = 99, seed=NULL,
                         order=1, mc=1,
                         test = c("rho", "cd", "sclm"),
                         alternative=c("twosided", "onesided",
                                       "symmetric"), \dots)
                          
}

\arguments{
  \item{x}{an object of class \code{formula}, \code{panelmodel}, or \code{pseries} 
          (depending on the respective interface) describing the model to be tested}
  \item{data}{a \code{data.frame}}
  \item{w}{a \code{n x n} \code{matrix} describing proximity between
    individuals, with \eqn{w_ij = a} where \eqn{a} is any number such
    that \code{as.logical(a)==TRUE},
    if \eqn{i,j} are neighbours, \eqn{0} or any number \eqn{b} such that
    \code{as.logical(b)==FALSE} elsewhere. Only the lower triangluar part
    (without diagonal) of \code{w} after coercing by \code{as.logical()}
    is evaluated for neighbouring information (but \code{w} can be
    symmetric). See also \bold{Details} and \bold{Examples}.}
 \item{index}{an optional numerical index, in case \code{data} has to be
  formatted by \code{plm.data}} 
 \item{model}{an optional character string indicating which type of
  model to estimate; 
   if left to \code{NULL}, the original heterogeneous specification of
   Pesaran is used}
 \item{replications}{the number of Monte Carlo randomizations of the
   neighbourhood matrix (default: 99),}
 \item{seed}{the optional random seed,}
 \item{order}{the order of neighbourhood to test for,}
 \item{mc}{the number of parallel threads to execute; defaults to 1
  (serial execution); is limited to the number of execution cores
  actually available, and depends on operating system support.}
 \item{test}{the type of test statistic to be returned. One of 
   \itemize{
    \item \code{"rho"} for the average correlation coefficient,
    \item \code{"cd"} for Pesaran's CD statistic, or
    \item \code{"sclm"} for the scaled version of Breusch and Pagan's LM
    statistic,}}
 \item{alternative}{the alternative hypothesis for the test, defaulting
  to (asymmetric) twosided,}
 \item{\dots}{further arguments to be passed on to \code{plm}, such as
  e.g. \code{effect} or \code{random.method}} 
}

\value{
An object of class \code{"htest"}.
}

\details{
  This test is meant as a generalization of Pesaran's spatial dependence
  test "CD(p)" for robustness against global dependence (perhaps of the
  factor type) and persistence in the data, both of which the original
  test does not tolerate. 

  The procedure can be applied to model residuals as well as to
  individual \code{pseries}.
  See the comments in \code{pcdtest} as for the different methods.
  
  Space is defined supplying a proximity matrix (elements coercible to
  \code{logical}) with argument \code{w} which provides information on
  whether any pair of individuals are neighbours or not. If
  \code{order=1}, only first-order neighbouring pairs will be used in
  computing the test; else, \code{w} will be transformed in the
  neighbourhood matrix of the appropriate order. The matrix need not be
  binary, so commonly used ``row--standardized'' matrices can be employed
  as well. \code{nb} objects from \pkg{spdep} must instead be transformed
  into matrices by \pkg{spdep}'s function \code{nb2mat} before using.

  Notice that the \code{"rho"} and \code{"cd"} tests are permutationally
  equivalent. 
 
  The test is suitable also for unbalanced panels.
 
  The test on a \code{pseries} is the same as a test on a pooled
  regression model of that variable on a constant,
  i.e. \code{rwtest(some_pseries)} is equivalent to
  \code{rwtest(plm(some_var ~ 1, data = some_pdata.frame, model =
  "pooling")} and also equivalent to \code{rwtest(some_var ~ 1, data =
  some_data)}, where \code{some_var} is the variable name in the data
  which corresponds to \code{some_pseries}.
}

\references{
 Millo, G. (2016), A simple randomization test for spatial dependence in
 the presence of common factors and serial correlation,
 \emph{(unpublished)}, \bold{xx}(x), pp. xxx--xxx. 
 Pesaran, M.H. (2004), General Diagnostic Tests for Cross Section
 Dependence in Panels, \emph{CESifo} Working Paper 1229. 
 Pesaran, M.H. (2015), Testing Weak Cross--Sectional Dependence in Large
 Panels, \emph{Econometric Reviews}, \bold{34}(6-10), pp. 1089--1117. 
}

\author{Giovanni Millo}

\examples{
data(Produc, package = "plm")
data(usaww)
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
## test on heterogeneous model (separate time series regressions)
rwtest(fm, data = Produc, w=usaww, index = c("state", "year"))

## test on two-way fixed effects homogeneous model
rwtest(fm, data = Produc, w=usaww, index = c("state", "year"),
       model = "within", effect = "twoways")

## test on panelmodel object
library(plm)
g <- plm(fm, data = Produc)
rwtest(g, w=usaww)

## test on pseries, higher-order neighbourhood
pprod <- pdata.frame(Produc)
rwtest(pprod$gsp, w=usaww, order=3)
}

\keyword{htest}
