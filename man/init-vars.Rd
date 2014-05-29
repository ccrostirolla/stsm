\name{init.vars}
\alias{init.vars}

\title{Initial Parameter Values}

\description{
This function computes initial variance parameters 
to be used as starting parameter values in an optimization procedure.
}

\usage{
init.vars(model, debug = FALSE)
}

\arguments{
\item{model}{an object of class \code{\link[stsm.class]{stsm}}.}
\item{debug}{logical. If \code{TRUE}, the correctness if the result is double-checked.}
}

\details{
%FIXME cite page or section of the book
As mentioned in Harvey (1989), the frequency domain representation of 
the structural model suggests using a linear regression to 
compute initial variance parameters from which to start an optimization procedure.
The variable \eqn{2\pi} times the periodogram is regressed on the constant terms
of the spectral generating function of the model.
}

\value{
A list containing the initial variance parameters and the output of
the linear regression.
}

\seealso{
\code{\link{maxlik.fd}},
\code{\link[stsm.class]{stsm.sgf}}.
}

\references{
Harvey, A. C. (1989).
\emph{Forecasting, Structural Time Series Models and the Kalman Filter}. 
Cambridge University Press.
}

\keyword{ts}
