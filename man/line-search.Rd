\name{linesearch}
\alias{linesearch}
\alias{Brent.fmin}
\alias{step.maxsize}

\title{Choice of the Step Size in the Scoring Algorithm}

\description{
These functions are used by the scoring optimization algorithms 
in order to choose the step size.
}

\usage{
Brent.fmin(a = 0, b, fcn, tol = .Machine$double.eps^0.25, ...)
linesearch(b, fcn, grd, ftol = 0.0001, gtol = 0.9, ...)
step.maxsize(x, xlo, xup, pd, cap = 1)
}

\arguments{
\item{a}{a numeric, lower end point of the interval where a line search procedure will 
search for the optimum step size. It should be zero or a positive value.}
\item{b}{a numeric, upper end point of the searching interval.}
\item{fcn}{a function to be optimized with respect to the step size.}
\item{tol}{the tolerance for convergence of the line search procedure.}
\item{grd}{a function returning the gradient of \code{fcn}.}
\item{ftol}{a numeric, tolerance of the Wolfe condition related to the value of the function.}
\item{gtol}{a numeric, tolerance of the Wolfe condition related to the gradient.}
\item{x}{a numeric containing the current value of the parameters.}
\item{xlo}{a numeric, lower bounds of the parameters of the model.}
\item{xup}{a numeric, upper bounds of the parameters of the model.}
\item{pd}{a numeric, direction vector chosen by the scoring algorithm.}
\item{cap}{the maximum step size allowed, the default is \eqn{1}.}
\item{...}{arguments to be passed to the objective function or the gradient.}
}

\details{
These functions are intended to be called by other functions, 
not to be used directly by the user. 
For details about how \code{fcn} and \code{grd} should be defined see 
the source code of \code{maxlik.fd.scoring} and \code{maxlik.td.scoring}.

The default line search procedure used by the scoring algorithms is the 
univariate optimization function \code{\link[stats]{optimize}} from 
package \code{stats}. 
The functions \code{linesearch} and \code{Brent.fmin} are used for debugging, 
didactic and experimental purposes. 
They provide useful information when testing 
the scoring algorithm and allowed easy tune of some parameters of the line search
procedure. This kind of information and options are not available for example 
in \code{\link[stats]{optim}} from package \pkg{stats}.

\code{Brent.fmin} is a version ported directly from the \R sources 
(procedure \code{Brent_fmin} in file \file{optimize.c}).
\code{linesearch} is based on Nocedal and Wright (2006) chapter 3 
and Pollock (1999) Chapter 12. It can be used 
to test the effect and role of the Wolfe conditions.

The function \code{step.maxsize} is not a line search procedure.
Given the direction vector chosen by the scoring algorithm,
this function returns the upper end of the interval where the line search procedure 
will search for the optimum step size. It ensures that for any step size inside the interval 
from \eqn{0} to the returned value the updated parameter values abide to the 
lower and upper bounds.
This approach is also used by A. Clausen in his implementation of the 
\href{http://economics.sas.upenn.edu/~clausen/computing/optim.php}{BFGS} algorithm.
The use of this function is a simple alternative to reparameterizations of the model 
and to the idea implemented in the L-BFGS-B algorithm in order to deal with this 
kind of constraints.
}

\value{
\code{Brent.fmin} and \code{linesearch} return a list containing:
\item{vx}{a vector containing the optimal value at each iteration during the bracketing.}
\item{minimum}{the optimal value found in the last iteration.}
\item{fx}{the value of the function for the optimal step size.}
\item{iter}{number of iterations employed by the procedure.}
\item{counts}{number of calls to the objective function. For \code{linesearch} it is 
a two-element vector where the second records the number of calls made to the 
gradient.}

\code{step.maxsize} returns a numeric containing the highest possible that is 
compatible with the arguments passed to the function (direction vector and bounds).
}

\seealso{
\code{\link[stsm]{maxlik.fd.scoring}},
\code{\link[stsm]{maxlik.td.scoring}}.
}

\references{
Brent, R. (1973) \emph{Algorithms for Minimization without Derivatives}.
Prentice-Hall.

Clausen, A. \R code for the BFGS algorithm.
\url{http://economics.sas.upenn.edu/~clausen/computing/optim.php}.

Nocedal, J. and Wright, J. W. (2006).
\emph{Numerical Optimization}. Springer-Verlag.

Pollock, D.S.G. (1999).
\emph{A Handbook of Time-Series Analysis Signal Processing and Dynamics}.
Academic Press.
}

\keyword{nonlinear}
\keyword{optimize}

