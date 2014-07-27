\name{maxlik.td}
\alias{maxlik.td}
\alias{maxlik.td.optim}
\alias{maxlik.td.scoring}

\title{Maximization of the Time Domain Likelihood Function}

\description{
Maximize the time domain log-likelihood function of a structural time series
model by means of a scoring algorithm or a general purpose optimization 
algorithm available in \code{\link[stats]{optim}}.
}

\usage{
maxlik.td.optim(m, 
  KF.version = eval(formals(KFKSDS::KalmanFilter)$KF.version),
  KF.args = list(), check.KF.args = TRUE,
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999, 
  method = c("BFGS", "L-BFGS-B", "Nelder-Mead", "CG", "SANN", "AB-NM"),
  gr = c("numerical", "analytical"), optim.control = list())

maxlik.td.scoring(m, step = NULL, 
  KF.args = list(), check.KF.args = TRUE,
  ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1),
  control = list(maxit = 100, tol = 0.001, trace = FALSE, silent = FALSE), 
  debug = FALSE)
}

\arguments{
\item{m}{an object of class \code{\link[stsm.class]{stsm}}.}
%
\item{KF.version}{character indicating the implementation of the Kalman filter to be used.}
%
\item{KF.args}{a list of parameters to be passed to the function chosen to run the Kalman filter.}
%
\item{check.KF.args}{logical. If \code{TRUE}, the arguments passed through \code{KF.args} are
checked for consistency with the interface chosen in \code{KF.version}:
for \code{maxlik.td.scoring} in it is checked for consistency with interface \code{KFKSDS}.}
%
\item{barrier}{a list defining a barrier term to penalize parameter values close to the bounds 
\code{m@lower} and \code{m@upper}.}
%
\item{inf}{a numeric indicating the value to be returned if the value of the log-likelihood 
function happens to be \code{NA} or non-finite at some iteration of the optimization procedure.}
%
\item{method}{character indicating the method to be used by \code{\link[stats]{optim}}.}
%
\item{gr}{character indicating whether numerical or analytical derivatives should be used.}
%
\item{optim.control}{a list of control parameters passed to \code{\link[stats]{optim}}.}
%
\item{step}{if it is a numeric it stands for a fixed step size, 
otherwise an automatic procedure is used to choose the step size.}
%
\item{ls}{control parameters for the line search procedure used to chose the step size
in the scoring algorithm.}
%
\item{control}{a list of control parameters for the scoring algorithm.}
%
\item{debug}{logical. If \code{TRUE}, tracing information is printed for debugging purposes 
of the scoring algorithm.}
}

\details{
The function \code{maxlik.td.optim} implements the common procedure of maximum likelihood, i.e.,
maximization of the time domain likelihood function by means of a numerical optimization 
algorithm (L-BFGS-B or other algorithms available in \code{\link[stats]{optim}}).
The likelihood function is evaluated by means of the Kalman filter.

The function \code{maxlik.td.scoring} implements a scoring algorithm based 
on the analytical expression of the information matrix of the time domain 
likelihood function.

López-de-Lacalle (2013) discusses several options and details that are often omitted 
or ignored when maximizing the likelihood function of a structural time series models.
The interface \code{maxlik.td.optim} allows the user to choose some options that 
may affect the results or convergence of the algorithm.

A novelty compared to other implementations such as \code{\link[stats]{StructTS}}
is that the optimization procedure is enhanced by means of analytical derivatives.
Another enhancement is that one of the parameters can be concentrated 
out of the likelihood function. The parameter to be concentrated is defined in 
the slot \code{cpar} of the input model \code{m}, 
see \code{\link[stsm.class]{stsm}} in package pkg{stsm.class}.
This option has not yet been implemented in the scoring algorithm.

For details about the options than can be passed through argument \code{KF.args}
see the documentation of the same argument in function \code{\link[KFKSDS]{KalmanFilter}}
in package \pkg{KFKSDS}.

For further information about the scoring algorithm see the following points in 
the details section of \code{\link{maxlik.fd.scoring}}:
\emph{Bounds on parameters and barrier term},
\emph{Control parameters for the scoring algorithm} and
\emph{Choice of the step size in the scoring algorithm}.

If external regressors are included in the model \code{m}, starting values for 
their coefficients are obtained in a linear regression of the differenced 
series on the differenced regressors. The values in the slot \code{pars} 
are therefore overwritten and not used as initial values.

Note: \code{ls$type = "wolfe"} is not implemented for \code{maxlik.td.scoring}.
}

\value{
A list of class \code{stsmFit}. 
See the section \sQuote{Value} in \code{\link{maxlik.fd}}.
}

\seealso{
\code{\link{barrier.eval}},
\code{\link{mloglik.td}},
\code{\link[stsm.class]{stsm-class}},
\code{\link[stats]{optim}}.
}

\references{
Durbin, J. and Koopman, S. J. (2001).
\emph{Time Series Analysis by State Space Methods}. 
Oxford University Press.

Harvey, A. C. (1989).
\emph{Forecasting, Structural Time Series Models and the Kalman Filter}. 
Cambridge University Press.

López-de-Lacalle, J. (2013).
\sQuote{101 Variations on a Maximum Likelihood Procedure for a Structural Time Series Model.}
Unpublished manuscript.

Nocedal, J. and Wright, J. W. (2006).
\emph{Numerical Optimization}. Springer-Verlag.
}

\examples{
# replicate maximum likelihood procedure as implemented in 'stats::StructTS'
res0 <- stats::StructTS(log(AirPassengers), type = "BSM")
mairp <- stsm.class::stsm.model(model = "BSM", y = log(AirPassengers), 
  transPars = "StructTS")
res1 <- maxlik.td.optim(mairp, KF.version = "KFKSDS", 
  KF.args = list(P0cov = TRUE), method = "L-BFGS-B", gr = "numerical")
mairp1 <- stsm.class::set.pars(mairp, pmax(res1$par, .Machine$double.eps))
round(stsm.class::get.pars(mairp1), 6)
all.equal(stsm.class::get.pars(mairp1), res0$coef[c(4,1:3)], 
  tol = 1e-04, check.attributes = FALSE)

# smoothed components
ss1 <- stsm.class::char2numeric(mairp1, P0cov = TRUE)
kf1 <- KFKSDS::KF(mairp1@y, ss1)
ks1 <- KFKSDS::KS(mairp1@y, ss1, kf1)
plot(tsSmooth(res0)[,c(1,3)])
plot(ks1$ahat[,c(1,3)])

# the scoring algorithm reaches another local optimum where
# the seasonal component is more homoscedastic
mairp <- stsm.class::stsm.model(model = "BSM", y = log(AirPassengers), 
  transPars = NULL)
res2 <- maxlik.td.scoring(mairp, step = NULL, 
  ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1),
  control = list(maxit = 100, tol = 0.001, trace = FALSE), debug = FALSE)
round(res2$par, 6)

mairp2 <- stsm.class::set.pars(mairp, res2$par)
ss2 <- stsm.class::char2numeric(mairp2, P0cov = FALSE)
kf2 <- KFKSDS::KF(mairp2@y, ss2)
ks2 <- KFKSDS::KS(mairp2@y, ss2, kf2)
plot(ks2$ahat[,c(1,3)])
}

\keyword{ts}
\keyword{models}
\keyword{nonlinear}
