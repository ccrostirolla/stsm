\name{maxlik.fd}
\alias{maxlik.fd}
\alias{maxlik.fd.optim}
\alias{maxlik.fd.scoring}
\alias{maxclik.fd.scoring}

\title{Maximization of the Spectral Likelihood Function}

\description{
Maximize the spectral log-likelihood function of a structural time series
model by means of a scoring algorithm or a general purpose optimization 
algorithm available in \code{\link[stats]{optim}}.
}

\usage{
maxlik.fd.optim(m, 
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999, 
  method = c("BFGS", "L-BFGS-B", "Nelder-Mead", "CG", "SANN"), 
  gr = c("analytical", "numerical"), optim.control = list())

maxlik.fd.scoring(m, step = NULL, 
  information = c("expected", "observed", "mix"),
  ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1),
  barrier = list(type = c("1", "2"), mu = 0), 
  control = list(maxit = 100, tol = 0.001, trace = FALSE, silent = FALSE), 
  debug = FALSE)
  
maxclik.fd.scoring(m, step = NULL, information = c("expected", "observed"),
  ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1),
  barrier = list(type = c("1", "2"), mu = 0), 
  control = list(maxit = 100, tol = 0.001, trace = FALSE, silent = FALSE))
}

\arguments{
\item{m}{an object of class \code{\link{stsm}}.}
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
\item{information}{the type of information about second order derivatives 
used to project the gradient in the scoring algorithm.}
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
The matrix used to project the gradient may be based on expected or observed 
information. The former, \code{information = "expected"}, uses the analytical expression of the 
information matrix. The latter, \code{information = "observed"}, uses the analytical
expression of the Hessian as in a Newton-Raphson procedure. The option \code{information = "mix"}
uses a mixture of both expressions; in simulations it performed similar to the 
information matrix, this option may be removed in future versions of the package.

\code{maxclik.fd.scoring} maximizes the concentrated likelihood function.
The parameter to be concentrated must be defined in the slot \code{cpar} of the input 
model \code{m}, see \code{\link{stsm}}.
\code{maxlik.fd.optim} detects whether \code{cpar} is defined in the input model.
In the scoring algorithm, if \code{m@cpar} is not \code{NULL} \code{maxclik.fd.scoring} 
should be used.

\emph{Bounds on parameters and barrier term.}
The lower and upper bounds within which the L-BFGS-B algorithm conducts the search 
are taken from the slots \code{lower} and \code{upper} defined in the input object \code{m}.
As an alternative to the L-BFGS-B procedure, a barrier term can be passed as argument.
%The argument \code{barrier} contains the following elements: 
%\code{type}, the type of barrier term used and \code{mu}, a positive value for the 
%barrier multiplier, or zero value if no barrier is considered.
The barrier term is added to the likelihood function and acts as a penalization 
for parameter values close to the bounds. 
For details about the barrier term see \code{\link[stsm]{barrier.eval}}.
In the scoring algorithm, if \code{step = NULL} the procedure automatically searches 
the optimum step size that is compatible with the bounds on the parameters. 
See \code{\link[stsm]{step.maxsize}} and further details below.

\emph{Control parameters for the scoring algorithm.}
\code{maxit}, maximum number of iterations (the default is \eqn{100});
\code{tol}, tolerance to assess convergence of the algorithm (the default is \eqn{0.001});
\code{trace}, logical, if \code{TRUE}, the values of the parameters at each iteration of 
the procedure are stored and returned in a matrix;
\code{silent}, logical, if \code{FALSE}, a warning is printed if convergence is not 
achieved.
%at the end of the process for the given tolerance and maximum number of iterations.

\emph{Choice of the step size in the scoring algorithm.}
If \code{step} is a numerical value, the step size is fixed to that value at all 
iterations of the algorithm. Otherwise, the choice of the step size in the scoring algorithm is 
made by means of a line search procedure specified in the argument \code{ls}.
See \code{\link[stsm]{step.maxsize}} for a description of the elements that can 
be passed through the argument \code{ls}.
%
%\itemize{
%\item{\code{type}, a character indicating the implementation of the line search procedure.
%The following are available: 
%\code{"optimize"}, which calls the function \code{\link[stats]{optimize}} from the \pkg{stats} package;
%\code{"brent.fmin"}, a version ported from \R sources (procedure \cpde{Brent_fmin} in file 
%\file{fmin.c}. This version is used for experimental and debugging purposes and may be removed 
%in future versions of the package;
%\code{"wolfe"}, based on Nocedal and Wright (2006) chapter 3.
%This version is also experimental and was used to check the effect of the Wolfe conditions.}
%%
%\item{\code{tol}, the tolerance for convergence of the line search procedure.}
%%
%\item{\code{cap}, the maximum step size allowed, the default is \eqn{1}.
%The actual maximum step size considered may be lower than the value \code{cap} passed in 
%\code{ls} if that step size leads to parameter values outside the lower and upper 
%bounds. In fact, choosing a step size at each step instead of defining it fixed 
%to the value \code{step} prevents the scoring algorithm from searching outside the 
%bounds. It is therefore an alternative to a parameterization for example in terms of 
%squares of the parameters that will ensure positive variances.}
%}

\emph{External regressors} 
If external regressors are included in the model \code{m}, starting values for 
their coefficients are obtained in a linear regression of the differenced 
series on the differenced regressors. The values in the slot \code{pars} 
are therefore overwritten and not used as initial values.
}

\value{
A list of class \code{stsmFit} with components:
\item{call}{an object of class \code{\link[base]{call}} specifying the arguments 
passed to the function.}
\item{init}{initial parameter values.}
\item{pars}{parameter values at the local optimum.}
\item{m}{the \code{\link{stsm}} model object 
updated with the optimal parameter values.}
\item{loglik}{the value of the log-likelihood function at the local optimum.}
\item{convergence}{convergence code returned by \code{\link[stats]{optim}}; for the 
scoring algorithm a logical indicating whether convergence was achieved.}
\item{iter}{for \code{maxlik.fd.optim} it is a two-element vector with the number of 
calls made by \code{\link[stats]{optim}} to the objective function and to the gradient; 
for \code{maxlik.fd.scoring} it is the number of iterations employed by the scoring algorithm.}
\item{message}{an empty character or a character message giving some additional information 
about the optimization process.}
\item{Mpars}{a matrix or \code{NULL}.
If \code{control$trace = TRUE} in the scoring algorithm, the path to the 
local optimum is traced storing by rows the parameter values at each iteration.}
\item{steps}{a vector or \code{NULL}. 
If \code{control$trace = TRUE} in the scoring algorithm, the step size used at each 
iteration is stored in this vector.}
\item{ls.iter}{a vector containing the number of iterations
employed by the line search procedure at each step of the scoring algorithm.
It is \code{NULL} if \code{ls$type = "optimize"}.}
\item{ls.counts}{a two-element vector containing the total number of calls to the 
objective function and the gradient made by the line search procedure in all the iterations.
(It is \code{NULL} if \code{ls$type = "optimize"}.)}

Note: if \code{m@transPars} is not \code{NULL}, the elements \code{init}
and \code{pars} are in terms of the auxiliary set of parameters. 
If the output is stored for example in an object called \code{res},
\code{get.pars(res$model)} will return the actual variance parameters.

The version based on \code{optim}, \code{maxlik.fd.optim}, returns also the element
\code{hessian} containing the numerically differentiated Hessian matrix 
at the local optimum.
Note that, if the model is parameterized in terms of an auxiliary set of parameters
and \code{gradient = "numerical"} is used, the Hessian returned by \code{optim} 
is defined with respect to the auxiliary set of parameters, not the variances. 
%FIXME see
%For the \sQuote{L-BFGS-B} method, the following note from \code{\link[stats]{optim}} applies,
%\emph{this is the Hessian of the unconstrained problem even if the box constraints are active}.
}

\seealso{
\code{\link{barrier.eval}},
\code{\link{mloglik.fd}},
\code{\link{stsm}},
\code{\link[stats]{optim}}.
}

\references{
Harvey, A. C. (1989).
\emph{Forecasting, Structural Time Series Models and the Kalman Filter}. 
Cambridge University Press.

Nocedal, J. and Wright, J. W. (2006).
\emph{Numerical Optimization}. Springer-Verlag.
}

\examples{
# fit the local level plus seasonal model to a 
# sample simulated series
# further examples and options can be explored in the 
# script files 'sim-llmseas-ml-fd.R' and 'sim-llmseas-mcl-fd.R'
# available in the 'inst' folder of the source package
data("llmseas")

# initial parameters and 'stsm' model
initpars <- c(var1 = 1, var2 = 1, var3 = 1)
m <- stsm.model(model = "llm+seas", y = llmseas, pars = initpars)

# Newton-Raphson algorithm (analytical Hessian)
res1 <- maxlik.fd.scoring(m = m, step = NULL, 
  information = "observed", control = list(maxit = 100, tol = 0.001))
res1

# Scoring algorithm (information matrix)
res2 <- maxlik.fd.scoring(m = m, step = NULL, 
  information = "expected", control = list(maxit = 100, tol = 0.001))
res2

# wrapper function for 'optim()' in the 'stats' package
res3 <- maxlik.fd.optim(m, method = "L-BFGS-B", gr = "analytical")
res3

# concentrating one of the parameters
# the model must be first defined accordingly, here 'var1', i.e.,
# the variance of the disturbance in the observation equation
# is concentrated, its standard error is reported as 'NA'
mc <- stsm.model(model = "llm+seas", y = llmseas, 
  pars = initpars[-1], cpar = initpars[1])
res4 <- maxclik.fd.scoring(m = mc, step = NULL, 
    information = "observed", control = list(maxit = 100, tol = 0.001))
res4
}

\keyword{ts}
\keyword{models}
\keyword{nonlinear}
