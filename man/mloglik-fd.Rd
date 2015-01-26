\name{mloglik.fd}
\alias{mloglik.fd}
\alias{mcloglik.fd}
\alias{mloglik.fd.deriv}
\alias{mcloglik.fd.deriv}
\alias{mloglik.fd.grad}
\alias{mcloglik.fd.grad}

\title{Spectral Log-Likelihood Function and Derivatives}

\description{
These functions evaluate the negative of the spectral log-likelihood function of 
a linear Gaussian state space model and its first and second order derivatives. 
}

\usage{
mloglik.fd(x, model, 
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999, xreg)

mcloglik.fd(x, model, xreg = NULL,
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999)

mloglik.fd.deriv(model, xreg = NULL,
  gradient = TRUE, hessian = TRUE, infomat = TRUE, modcovgrad = TRUE,
  barrier = list(type = c("1", "2"), mu = 0),
  version = c("2", "1"))

mcloglik.fd.deriv(model, xreg = NULL, 
  gradient = TRUE, hessian = TRUE, infomat = TRUE)

mloglik.fd.grad(x, model, xreg = NULL,
  barrier = list(type = c("1", "2"), mu = 0),
  inf)

mcloglik.fd.grad(x, model, xreg = NULL, inf, barrier)
}

\arguments{
\item{x}{a numeric vector containing the parameters of the model.
This is an auxiliar argument so that this function can be used as input
to \code{\link[stats]{optim}}.}
%
\item{model}{an object of class \code{\link{stsm}}.}
%
\item{xreg}{optional list containing constant terms. See details.}
%
\item{barrier}{a list defining a barrier term to penalize parameter values close to the bounds 
\code{m@lower} and \code{m@upper}.}
%
\item{inf}{a numeric indicating the value to be returned if the value of the log-likelihood 
function happens to be \code{NA} or non-finite.}
%
\item{gradient}{logical. If \code{TRUE}, first order derivatives of the
negative of the spectral log-likelihood function are returned.}
%
\item{hessian}{logical. If \code{TRUE}, second order derivatives of the
negative of the spectral log-likelihood function are returned.}
%
\item{infomat}{logical. If \code{TRUE}, the information matrix
of the spectral log-likelihood are returned.}
%
\item{modcovgrad}{logical. If \code{TRUE}, a mixture of the analytical expressions 
for the Hessian and the outer product of the gradient is used.
This option is experimental and may be removed in future versions of the package.
}
%
\item{version}{a character indicating whether implementation \code{"2"} or \code{"2"}
(the default) should be used.They yield the same result but are kept for debugging and 
comparison of timings. This argument may be removed in future versions.}
}

\details{
The spectral log-likelihood of a linear Gaussian state space model is given by
(Harvery, 1989 Section 4.3):

\deqn{logLik = -0.5 log(2\pi) - 0.5 \sum_{j=0}^{n-1} log\,g(\lambda[j]) -
\pi \sum_{j=0}^{n-1} I(\lambda[j])/g(\lambda[j])}

where \eqn{\lambda[j]} is a frequency defined as
\eqn{\lambda[j] = 2\pi j/n}; 
\eqn{I(\lambda[j])} is the periodogram at frequency \eqn{\lambda[j]}
and
\eqn{g(\lambda[j])} is the spectral generating function of the model
at frequency \eqn{\lambda[j]}.

The derivation of the spectral likelihood function defined above relies 
on the assumption that the process is circular (its covariance matrix 
is circulant). If the process is not circular the value of
the likelihood is an approximation.

First and second order derivatives are computed by means of their 
analytical expressions.
The first order derivatives of the spectral log-likelihood with respect to 
parameter \eqn{\theta} are given by:
\deqn{
(d\,logLik)/(d\,\theta) =
0.5 \sum_{j=0}^{n-1} ( (2\pi I(\lambda[j])) / (g(\lambda[j])) - 1 )
(1/g(\lambda[j])) d g(\lambda[j])/d \theta}

Second order derivatives are given by:
\deqn{
(d^2\, logLik)(d\,\theta \theta') =
\sum_{j=0}^{n-1} ( (2\pi I(\lambda[j]))/(g(\lambda[j])) - 1 ) 
1/(2g(\lambda[j])) (d^2 g(\lambda[j]))/(d \theta d \theta') -}
\deqn{
2 \sum_{j=0}^{n-1} ( (4\pi I(\lambda[j]))(g(\lambda[j])) - 1 )
( \frac{1}{2g(\lambda[j])} )^2 
(d g(\lambda[j]))/(d \theta)
(d g(\lambda[j]))/(d \theta')
}

The argument \code{x} is an auxiliar vector that is necessary in some contexts. 
For example, the input to function \code{\link[stats]{optim}} must contain
as first argument the vector of parameters where optimization is performed.
If it is not required or is redundant information contained in \code{model@pars} 
it can be set to \code{NULL}.

The functions \code{mcloglik.fd}, \code{mcloglik.fd.deriv} and 
\code{mcloglik.fd.grad} use the expressions for the spectral log-likelihood 
function where the parameter specified in \code{model@cpar} is concentrated
out of the likelihood function.

For further information about the barrier term see 
\emph{Bounds on parameters and barrier term}
in the details section in \code{\link{maxlik.fd.scoring}}.

Arguments \code{inf} and \code{barrier} are not used by \code{mloglik.fd.grad} 
and \code{mcloglik.fd.grad} but they are needed in \code{\link{maxlik.fd.optim}}, 
where this function is passed as the gradient to be used by \code{\link[stats]{optim}}
including the arguments \code{inf} and \code{barrier}.
%\code{inf} and \code{barrier} are in turn passed as argument 
%to the objective function \code{\link{mloglik.fd}}.

Argument \code{xreg}.
It is an optional list of constant terms. It is used by 
\code{\link{maxlik.fd.optim}} when analytical derivatives are employed and 
by \code{\link{maxlik.fd.scoring}}. It avoids computing some constant terms each 
time the function \code{mloglik.fd.grad} is called. 

The list \code{xreg} should contain 
an element called \code{dxreg}, the external regressors differeneced by means 
of the differencing filter that renders stationarity in the model and 
the element \code{fft.xreg}, the Fourier transform of each regressor 
in \code{dxreg}.

The list \code{xreg} is not used by \code{mloglik.fd}. It is necessary 
to define this argument in the prototype of the function because when this function 
is passed to \code{\link[stats]{optim}} along with 
\code{mloglik.fd.grad}, the argument \code{xreg} is passed to \code{mloglik.fd}
when it is defined in \code{optim} as an argument to be passed 
to \code{mloglik.fd.grad}.

Argument \code{xreg} is not currently implemented in functions with 
concentration of a parameter, \code{mloglik.fd}, \code{mcloglik.fd.deriv}
and \code{mcloglik.fd.grad}.

Note: \code{modcovgrad} is not available when external regressors are defined 
in the input model, \code{model}.
}

\note{
\code{mcloglik.fd.deriv} is not currently implemented for model with non-null 
\code{model@transPars} or with a barrier term.
}

\value{
\code{mloglik.fd} returns a numeric value of the negative 
of the spectral log-likelihood evaluated at the parameter values defined in 
the input model \code{model} or at \code{x} if this argument is not \code{NULL}.
If the value happens to be \code{NA} or non-finite the value of argument \code{inf} 
is returned. This function is suited to be passed as the objective function to 
\code{\link[stats]{optim}}.

\code{mloglik.fd.deriv} returns a list containing a vector of the first order derivatives 
of the negative of the spectral likelihood and a matrix for the second order 
derivatives. Those derivative terms that are not requested by setting the 
corresponding argument to \code{FALSE} are set to \code{NULL} in the output list.

\code{mloglik.fd.grad} returns a numeric vector containing the gradient.
This function is suited to be passed as the gradient function to 
\code{\link[stats]{optim}}.

\code{mcloglik.fd}, \code{mcloglik.fd.deriv} and \code{mcloglik.fd.grad}
return the value of the same information as the other functions but for the 
concentrated likelihood function.
}

\seealso{
\code{\link{barrier.eval}},
\code{\link{logLik}},
\code{\link{maxlik.fd}},
\code{\link{stsm}}.
}

\references{
Harvey, A. C. (1989).
\emph{Forecasting, Structural Time Series Models and the Kalman Filter}.
Cambridge University Press.
}

\examples{
# define the local level model for Nile time series
pars <- c("var1" = 11000, "var2" = 1700)
m <- stsm.model(model = "local-level", y = Nile, pars = pars)
# 'mloglik.fd' returns the negative of the log-likelihood function
mloglik.fd(model = m)
# 'logLik' returns the value of the log-likelihood function
logLik(object = m, domain = "frequency")

# compare analytical and numerical derivatives
# more tests in file 'test-derivatives-mloglik-fd.R' in the 
# folder 'inst' of the source package
system.time(da <- mloglik.fd.deriv(m, gradient = TRUE, hessian = TRUE))
dgn <- numDeriv::grad(func = mloglik.fd, x = m@pars, model = m)
dhn <- numDeriv::hessian(func = mloglik.fd, x = m@pars, model = m)
all.equal(as.vector(da$gradient), dgn)
all.equal(da$hessian, dhn)

# the same as above for the local level plus seasonal model and 
# a sample simulated series
data("llmseas")
m <- stsm.model(model = "llm+seas", y = llmseas)
system.time(a <- mloglik.fd.deriv(model = m, gradient = TRUE, hessian = TRUE))
system.time(g <- numDeriv::grad(func = mloglik.fd, x = m@pars, model = m))
system.time(h <- numDeriv::hessian(func = mloglik.fd, x = m@pars, model = m))
all.equal(a$gradient, g, check.attributes = FALSE)
all.equal(a$hessian, h, check.attributes = FALSE)
}

\keyword{ts}
