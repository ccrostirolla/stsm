\name{barrier.eval}
\alias{barrier.eval}

\title{Barrier Term in the Likelihood Function}

\description{
This function evaluates a barrier term to be added to 
the log-likelihood function in order to penalize for parameter 
values close to the boundaries.
}

\usage{
barrier.eval(model, type = c("1", "2"), mu = 0.01, 
  gradient = FALSE, hessian = FALSE)
}

\arguments{
\item{model}{object of class \code{\link[stsm.class]{stsm}}.}
\item{type}{a character indicating the type of barrier term.}
\item{mu}{the barrier multiplier. A positive value.}
\item{gradient}{logical. If \code{TRUE}, the first order derivatives of the
barrier term with respect to the parameters of the model are evaluated.}
\item{hessian}{logical. If \code{TRUE}, the second order derivatives of the
barrier term with respect to the parameters of the model are evaluated.}
}

\details{
%Two types of \bold{barrier terms} are considered in order to penalize against 
%parameter values beyond the boundaries given in \code{lower} and \code{upper}.
%\deqn{-log(pars - lower) or 1/(pars - lower)}
%for constraints of the type \eqn{pars \ge lower}; and
%\deqn{-log(upper - pars) or 1/(upper - pars)}
%for constraints of the type \eqn{pars \le lower}.
%The logarithm barrier term is added in \code{FD.mloglik.barrier1} and 
%the latter barrier term is added by \code{FD.mloglik.barrier2}.
%
%The inequality constraints given in \code{model@lower} and \code{model@upper} 
%are dropped by means of a barrier term that is added to the 
%minus log-likelihood function.

A barrier term can be defined in order to penalize against 
parameter values close to the boundaries defined in \code{model@lower} 
and \code{model@upper}. The barrier term is added to the 
negative of the log-likelihood function, which is then given by:

\deqn{-logLik + \mu \sum_i q_i(x).}

\eqn{\mu} is a positive barrier multiplier, the larger it is the stronger 
the barrier is and, hence, the search is done farther from the boundaries
(Rardin, 1998 Fig. 14.8);
\eqn{q_i(x)} are the barrier functions that are defined in such a way that 
the barrier term grows to infinity as the parameters of the model approach to
the boundaries. Two types of barrier functions are considered:

\describe{
\item{\code{type = "1"}:}{\eqn{-log(par[i] - bound[i])}, for lower bound constraints.}
\item{\code{type = "1"}:}{\eqn{-log(bound[i] - par[i])}, for upper bound constraints.}
\item{\code{type = "2"}:}{\eqn{1/(par[i] - bound[i])}, for lower bound constraints.}
\item{\code{type = "2"}:}{\eqn{1/(bound[i] - par[i])}, for upper bound constraints.}
}
}

\value{
A list containing the barrier term evaluated for the lower and upper bound 
constraints, \code{bl} and \code{bu}, respectively, and the first and second 
order derivatives of the barrier term, \code{dl1}, \code{du1}, 
\code{dl2}, \code{du2}.
}

\seealso{
\code{\link{mloglik.fd}},
\code{\link{mloglik.td}}.
}

\references{
Rardin, R.L. (1998). Section 14.5. 
\emph{Optimization in Operations Research}. Prentice Hall.
}

\examples{
# define local level plus seasonal model for a simulated time series
# and evaluate the barrier term for some parameter values
# the barrier term in the second model is higher since the variances 
# are closer to the lower bounds (zero)
data("llmseas")
pars <- c("var1" = 30, "var2" = 1, "var3" = 10)
m1 <- stsm.model(model = "llm+seas", y = llmseas[,1], pars = pars)
bar1 <- barrier.eval(m1, type = "2", mu = 3)
bar1$barrier

pars <- c("var1" = 0.3, "var2" = 0.01, "var3" = 0.1)
m2 <- stsm.model(model = "llm+seas", y = llmseas[,1], pars = pars)
bar2 <- barrier.eval(m2, type = "2", mu = 3)
bar2$barrier

# the barrier term is added to the negative of the likelihood function
# that will be the objective function to be minimized,
# value of minus log-likelihood without barrier
mloglik.fd(model = m2)
# adding a barrier term
mloglik.fd(model = m2) + bar2$barrier
mloglik.fd(model = m2, barrier = list(type = "2", mu = 3))
}

\keyword{nonlinear}
