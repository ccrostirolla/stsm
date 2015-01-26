\name{stsm.sgf}
\alias{stsm.sgf}

\title{Spectral Generating Function of Common Structural Time Series Models}

\description{
Evaluate the spectral generating function of  of common
structural models: local level model, local trend model and
basic structural model.
}

\usage{
stsm.sgf(x, gradient = FALSE, hessian = FALSE, deriv.transPars = FALSE)
}

\arguments{
\item{x}{object of class \code{\link{stsm}}.}

\item{gradient}{logical. If \code{TRUE}, the gradient is returned.}
\item{hessian}{logical. If \code{TRUE}, \code{hessian} the gradient is returned.}
\item{deriv.transPars}{logical. If \code{TRUE}, the gradient and the Hessian are 
scaled by the gradient of the function that transforms the parameters.
Ignored if \code{x@transPars} is null.}
}

\details{
The stationary form of the \bold{local level model} is 
(\eqn{\Delta} is the differencing operator):

\deqn{\Delta y[t] = v[t] + \Delta e[t]}

and its \emph{spectral generating function} 
at each frequency \eqn{\lambda[j] = 2\pi j/T} for \eqn{j=0,...,T-1} is:

\deqn{g(\lambda[j]) = \sigma^2_2 + 2(1 - cos \lambda[j]) \sigma^2_1}

The stationary form of the \bold{local trend model} for a time series of
frequency \eqn{S} is:

\deqn{\Delta^2 y[t] = \Delta v[t] + w[t-1] + \Delta^2 e[t]}

and its \emph{spectral generating function} is:

\deqn{g(\lambda[j]) = 2(1 - cos \lambda[j]) \sigma^2_2 + \sigma^2_3 + 
4(1 - cos \lambda[j]) \sigma^2_1}

The stationary form of the \bold{basic structural model} for a time series
of frequency \eqn{p} is:

\deqn{\Delta \Delta^p y[t] = \Delta^p v[t] + S(L) w[t-1] + 
\Delta^2 s[t] + \Delta \Delta^p e[t]}

and its \emph{spectral generating function} is:

\deqn{g(\lambda[j]) = g_v(\lambda[j]) \sigma^2_2 + g_w(\lambda[j]) \sigma^2_3 + 
g_s(\lambda[j]) \sigma^2_4 + g_e(\lambda[j]) \sigma^2_1}

with

\deqn{g_v(\lambda[j]) = 2(1 - cos(\lambda[j] p))}
\deqn{g_w(\lambda[j]) = (1 - cos(\lambda[j] p))/(1 - cos(\lambda[j]))}
\deqn{g_s(\lambda[j]) = 4 (1 - cos(\lambda[j]))^2}
\deqn{g_e(\lambda[j]) = 4 (1 - cos(\lambda[j])) (1 - cos(\lambda[j] p))}
}

\value{
A list containing the following results:

\item{sgf}{spectral generating function of the BSM model at
each frequency \eqn{\lambda[j]} for \eqn{j=0,\dots,T-1}.}
\item{gradient}{first order derivatives of the spectral generating function with
respect of the parameters of the model.}
\item{hessian}{second order derivatives of the spsectral generating function with
respect of the parameters of the model.}
\item{constants}{the terms \eqn{g_v(\lambda[j])}, \eqn{g_w(\lambda[j])},
\eqn{g_s(\lambda[j])} and \eqn{g_e(\lambda[j])} 
that do not depend on the variance parameters.}
}

\seealso{
\code{\link{set.sgfc}},
\code{\link{stsm}},
\code{\link{stsm.model}}.
}

\references{
Harvey, A. C. (1989).
\emph{Forecasting, Structural Time Series Models and the Kalman Filter}. 
Cambridge University Press.
}

\examples{
# spectral generating function of the local level plus seasonal model
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 15), nopars = c("var3" = 30))
res <- stsm.sgf(m)
res$sgf
plot(res$sgf)
res$constants
# the element 'constants' contains the constant variables
# related to each component regardless of whether the 
# variances related to them are in the slot 'pars' or 'nopars'
names(get.pars(m))
colnames(res$constants)

# compare analytical and numerical derivatives
# identical values
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 15, "var3" = 30))
res <- stsm.sgf(m, gradient = TRUE)

fcn <- function(x, model = m) {
  m <- set.pars(model, x)
  res <- stsm.sgf(m)
  sum(res$sgf)
}

a1 <- numDeriv::grad(func = fcn, x = get.pars(m))
a2 <- colSums(res$grad)
all.equal(a1, a2, check.attributes = FALSE)

# analytical derivatives are evaluated faster than numerically
system.time(a1 <- numDeriv::grad(func = fcn, x = get.pars(m)))
system.time(a2 <- colSums(stsm.sgf(m, gradient = TRUE)$grad))
}

\keyword{ts}
