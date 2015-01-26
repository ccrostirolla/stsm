\name{method-logLik}
\alias{method-logLik}
\alias{logLik}
\alias{logLik.stsm}

\title{Extract Log-Likelihood}

\description{
This method returns the log-likelihood value of a model defined 
in a \code{\link{stsm}} object.
}

\usage{
\method{logLik}{stsm}(object, domain = c("frequency", "time"), xreg = NULL,
  td.args = list(P0cov = FALSE, t0 = 1, 
    KF.version = eval(formals(KFKSDS::KalmanFilter)$KF.version)), 
  check.td.args = TRUE, 
  barrier = list(type = c("1", "2"), mu = 0), 
  inf = 99999, ...)
}

\arguments{
\item{object}{an object of class \code{\link{stsm}}.}
\item{domain}{a character indicating whether the time domain or the spectral 
likelihood should be evaluated.}
\item{xreg}{optional matrix of external regressors with the same number of rows 
as the length of the input time series \code{m@y}.}
\item{td.args}{a list containing arguments to be passed to \code{\link{mloglik.td}}.
Only for \code{domain = "time"}.}
\item{check.td.args}{logical. If \code{TRUE}, argument \code{td.args} is checked 
for consistency with \code{td.args$KF.version}. Only for \code{domain = "time"}.}
\item{barrier}{a list defining a barrier term to penalize parameter values close to the bounds 
\code{m@lower} and \code{m@upper}.}
\item{inf}{a numeric indicating the value to be returned if the value of the log-likelihood 
function happens to be \code{NA} or non-finite.}
\item{...}{additional arguments. Currently ignored.}
}

\value{
An object of class \code{logLik} containing the value of the log-likelihood function 
for the given model and the attribute \code{df}, the number of estimated parameters.
}

\seealso{
\code{\link{mloglik.fd}},
\code{\link{mloglik.td}},
\code{\link{mle}},
\code{\link[KFKSDS]{KalmanFilter}},
\code{\link{stsm}}.
}

\keyword{ts}
