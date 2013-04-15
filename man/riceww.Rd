\name{riceww}
\docType{data}
\alias{riceww}
\title{Spatial weights matrix of Indonesian rice farms}

\description{
Spatial weights matrix of the 171 farms in the Indonesian Rice Farming example. Farms in the same village (out of six) are considered contiguous.  

}

\usage{
data(riceww)
}

\format{
A matrix with elements different from zero if farms i and j are
neighbors. Farms are considered neighbors if in the same village. Weights are row-standardized.}

\author{Giovanni Millo, data provided by Yves Croissant}

\keyword{datasets}
