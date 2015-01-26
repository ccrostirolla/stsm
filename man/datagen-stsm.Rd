\name{datagen.stsm}
\alias{datagen.stsm}

\title{Generate Data from a Structural Time Series Model}

\description{
This function simulates data from a structural time series model 
defined in an object of class \code{\link{stsm}}.
}

\usage{
datagen.stsm(n, model = list(), SigmaEV, labels, n0 = 20, freq = 1, 
  old.version = FALSE)
}

\arguments{
\item{n}{number of observations in the output time series.}
\item{model}{a list containing the matrices of the state space for of the 
structural model.}
\item{SigmaEV}{a list containing the elements \code{values} and \code{vectors},
they are respectively the eigen values and vectors of the covariance matrix of the 
disturbance terms in the state equation of the model.}
\item{labels}{optional vector of characters giving the names of the 
unobserved components in the model.}
\item{n0}{number of warming-up observations (they are not included in the output data).}
\item{freq}{number of observations per unit of time, e.g. \code{freq = 4} for 
quarterly data.}
\item{old.version}{logical. If \code{TRUE}, results obtained in a previous version 
of the package are reproduced.}
}

\details{
The matrices in the argument \code{model} must follow the conventions of an object
of class \code{\link{stsm}} as returned by
\code{\link{char2numeric}}.

For compatibility with previous versions of the package, 
\code{old.version = TRUE} generates random values from the multivariate normal distribution
using the function \code{\link[mvtnorm]{rmvnorm}} with \code{pre0.9_9994 = TRUE}.
\code{mvrnorm.version = FALSE} uses the theoretical expression that is commonly 
used to define random draws for the multivariate normal distribution.
}

\value{
A list containing the output time series and the underlying components 
of the model.
}

\seealso{
\code{\link{sim-data}},
\code{\link{stsm}}.
}

\examples{
# generate a quarterly series from a local level plus seasonal model
# the data set 'llmseas' is generated as follows (first series)
pars <- c(var1 = 300, var2 = 10, var3 = 100)
m <- stsm.model(model = "llm+seas", y = ts(seq(120), frequency = 4), 
  pars = pars, nopars = NULL)
ss <- char2numeric(m)
set.seed(123)
y <- datagen.stsm(n = 120, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q), 
  n0 = 20, freq = 4, old.version = TRUE)$data

data("llmseas")
all.equal(y, llmseas)
}

\keyword{ts}
