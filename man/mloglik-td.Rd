\name{mloglik.td}
\alias{mloglik.td}
\alias{KFconvar}
\alias{mloglik.td.deriv}
\alias{mloglik.td.grad}

\title{Time Domain Log-Likelihood Function and Derivatives}

\description{
This function evaluates the negative of the time domain log-likelihood function 
of a linear Gaussian state space model by means of the Kalman filter.
}

\usage{
mloglik.td(x, model, xreg = NULL,
  KF.version = eval(formals(KFKSDS::KalmanFilter)$KF.version),
  KF.args = list(), check.KF.args = TRUE,
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999)

KFconvar(model, P0cov = FALSE, barrier = list(type = "1", mu = 0), debug = TRUE)

mloglik.td.deriv(model, gradient = TRUE, infomat = TRUE,
  KF.args = list(), version = c("1", "2"), kfres = NULL)

mloglik.td.grad(x, model, KF.version, KF.args = list(), check.KF.args, barrier, inf, xreg)
}

\arguments{
\item{x}{a numeric vector containing the parameters of the model.
This is an auxiliar argument so that this function can be used as input
to \code{\link[stats]{optim}}.}
%
\item{model}{an object of class \code{\link[stsm.class]{stsm}}.}
%
\item{xreg}{optional matrix of external regressors with the same number of rows 
as the length of the input time series \code{m@y}.}
%
\item{KF.version}{character indicating the implementation of the Kalman filter to be used.}
%
\item{KF.args}{a list of parameters to be passed to the function chosen to run the Kalman filter.}
%
\item{check.KF.args}{logical. If \code{TRUE}, the elements passed in argument \code{KF.args} are 
checked for consistency with \code{KF.version}.}
%
\item{barrier}{a list defining a barrier term to penalize parameter values close to the bounds 
\code{m@lower} and \code{m@upper}.}
%
\item{inf}{a numeric indicating the value to be returned if the value of the log-likelihood 
function happens to be \code{NA} or non-finite.}
%
\item{P0cov}{logical. If \code{TRUE}, values outside the diagonal of the covariance matrix 
of the initial state vector are set equal to the values in the diagonal,
as done in \code{\link[stats]{StructTS}}}.
%
\item{debug}{logical. If \code{TRUE}, the correctness of the result is double-checked.}
%
\item{gradient}{logical. If \code{TRUE}, first order derivatives of the
negative of the spectral log-likelihood function are returned.}
%
\item{infomat}{logical. If \code{TRUE}, the information matrix
of the spectral log-likelihood are returned.}
%
\item{version}{a character indicating whether implementation \code{"2"} or \code{"2"}
(the default) should be used.They yield the same result but are kept for debugging and 
comparison of timings. This argument may be removed in future versions.}
%
\item{kfres}{optional list containing the elements involved in the Kalman filter 
as returned by \code{\link[KFKSDS]{KF.deriv}}.}
}

\details{
The general univariate linear Gaussian state space model is defined 
as follows:
\deqn{y[t] = Za[t] + e[t], e[t] \sim N(0, H)}
\deqn{a[t+1] = Ta[t] + Rw[t], w[t] \sim N(0, V)}

for \eqn{t=1,\dots,n} and \eqn{a[1] \sim N(a0, P0)}.
\eqn{Z} is a matrix of dimension \eqn{1\times m}{1xm}; 
\eqn{H} is \eqn{1\times 1}{1x1};
\eqn{T} is \eqn{m\times m}{mxm};
\eqn{R} is \eqn{m\times r}{mxr};
\eqn{V} is \eqn{r\times r}{rxr}; 
\eqn{a0} is \eqn{m\times 1}{mx1} and
\eqn{P0} is \eqn{m\times m}{mxm},
where \eqn{r} is the number of variance parameters in the state vector.

The Kalman filtering recursions for the model above are:

\emph{Prediction}
\deqn{a[t] = T a[t-1]}
\deqn{P[t] = T P[t-1] T' + R V R'} 
\deqn{v[t] = y[t] - Z a[t]}
\deqn{F[t] = Z P[t] Z' + H} 

\emph{Updating}
\deqn{K[t] = P[t] Z' F[t]^{-1}}
\deqn{a[t] = a[t] + K[t] v[t] }
\deqn{P[t] = P[t] - K[t] Z P[t]'}

for \eqn{t=2,\dots,n}, starting with \eqn{a[1]} and \eqn{P[1]} equal
to \code{a0} and \code{P0}. \eqn{v[t]} is the prediction error at observation 
in time \eqn{t} and \eqn{F[t]} is the variance of \eqn{v[t]}.

The log-likelihood of the model for a given set of parameter values is:

\deqn{logLik = -0.5 log(2\pi) - 0.5 \sum_{t=1}^n { log F[t] + v[t]^2 / F[t] }}

For details about the options than can be passed through argument \code{KF.args}
see the documentation of the same argument of function 
\code{\link[KFKSDS]{KalmanFilter}} in package \pkg{KFKSDS}.
For \code{mloglik.td.deriv}, the only element that is used if provided in 
\code{KF.args} is \code{P0cov} (a logical indicating whether the covariance matrix 
of the initial state vector diagonal or not).

The argument \code{x} is an auxiliar vector that is necessary in some contexts. 
For example, the input to function \code{\link[stats]{optim}} must contain
as first argument the vector of parameters where optimization is performed.
If it is not required or is redundant information contained in \code{model@pars} 
it can be set to \code{NULL}.

For further information about the barrier term see 
\emph{Bounds on parameters and barrier term}
in the details section in \code{\link[stsm]{maxlik.fd.scoring}}.

\code{KFconvar} evaluates the concentrated likelihood function.
The likelihood is concentrated with respect to the parameter defined in \code{model@cpar}.
The optimal value of the parameter that is concentrated out of the likelihood is:
%
\deqn{s2 = (1/n) \sum_{t=1}^n v[t]/ F[t]}
%
and the concentrated likelihood function is given by:
%
\deqn{clogLik = (n/2) log(2\pi + 1) + 0.5 \sum_{t=1}^n log(f[t]) + (n/2) log(s2).}

The gradient and the information matrix are calculated upon their corresponding
analytical expressions.

Arguments \code{KF.version}, \code{check.KF.args}, \code{barrier} and \code{inf} 
are not used by \code{mloglik.td.grad} 
but they are needed in \code{\link{maxlik.td.optim}}, where this function is passed 
as the gradient to be used by \code{\link[stats]{optim}} 
including the arguments \code{inf} and \code{barrier}.
%and these arguments are in turn passed 
%to the objective function \code{\link{mloglik.td}}.
}

\value{
The minus log-likelihood function evaluated at the parameter values defined in 
the input model \code{model} or at \code{x} if this argument is not \code{NULL}.
If the value happens to be \code{NA} or non-finite the value of argument 
\code{inf} is returned.
This function is suited to be passed as the objective function to 
\code{\link[stats]{optim}}.

\code{KFconvar} returns a list containing the element \code{mll}, 
the negative of the concentrated minus log-likelihood function and
the element \code{cpar}, the optimal value of the parameter that is 
concentrated out of the likelihood.

\code{mloglik.td.deriv} returns a list containing a vector of the first order derivatives 
of the negative of the time domain likelihood function and a matrix for the information 
matrix. They are set to \code{NULL} if any of them are not requested.

\code{mloglik.td.grad} returns a numeric vector containing the gradient.
This function is suited to be passed as the gradient function to 
\code{\link[stats]{optim}}.
}

\seealso{
\code{\link{barrier.eval}},
\code{\link{logLik}},
\code{\link{maxlik.td}},
\code{\link[stsm.class]{stsm-class}},
\code{\link[KFKSDS]{KalmanFilter}},
\code{\link[stats]{KalmanLike}}.
}

\references{
Durbin, J. and Koopman, S. J. (2001).
\emph{Time Series Analysis by State Space Methods}. 
Oxford University Press.

Harvey, A. C. (1989).
\emph{Forecasting, Structural Time Series Models and the Kalman Filter}. 
Cambridge University Press.
}

\examples{
# local level plus seasonal model for a sample simulated series
data("llmseas")
m <- stsm.class::stsm.model(model = "llm+seas", y = llmseas,
  pars = c("var1" = 300, "var2" = 10, "var3" = 100))
# evaluate the time domain likelihood function using 
# excluding the contributions of the first 8 observations
mloglik.td(model = m, KF.version = "KFKSDS", KF.args = list(t0 = 9))

# compare analytical and numerical derivatives
# identical gradient up to a tolerance
a <- mloglik.td.deriv(m, infomat = TRUE)
g <- numDeriv::grad(func = mloglik.td, x = m@pars, 
  model = m, KF.version = "KFKSDS")
h <- numDeriv::hessian(func = mloglik.td, x = m@pars, 
  model = m, KF.version = "KFKSDS")
all.equal(a$gradient, g, check.attributes = FALSE)
}

\keyword{ts}
