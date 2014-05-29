\name{maxlik.em}
\alias{maxlik.em}

\title{Maximization of the Time Domain Likelihood Function 
via the Expectation-Maximization Algorithm}

\description{
Maximize the time domain log-likelihood function of a structural time series
model by means of the Expectation-Maximization algorithm.
}

\usage{
maxlik.em(model, type = c("standard", "modified", "mix"), 
  tol = 0.001, maxiter = 300, kfconv = c(0, 10, 1),  
  ur.maxiter = 1000, r.interval = c(0.001, var(model@y)*0.75, 20),
  mod.steps = seq(3, maxiter, 10), 
  parallel = FALSE, num.cores = NULL)
}

\arguments{
\item{model}{an object of class \code{\link[stsm.class]{stsm}}.}
\item{type}{a character choosing the implementation of the algorithm, 
the original (standard), a modified version or a mixture of both.}
\item{tol}{the convergence tolerance.}
\item{maxiter}{the maximum number of iterations.}
\item{kfconv}{a length-three vector with control parameters used 
to determine whether the Kalman filter has converged at a given iteration.}
\item{ur.maxiter}{the maximum number of iterations used in the root finding procedure.
Ignored with \code{type = "standard"}.}
\item{r.interval}{a three-length vector to control the ends of the interval in the 
root finding procedure.}
\item{mod.steps}{a vector indicating the iterations at which the modified EM algorithm 
is to be run. Only for \code{type = "mix"}.}
\item{parallel}{logical. If \code{TRUE} the process is run in parallel by means of function 
\code{\link[parallel]{mclapply}} in package \pkg{parallel}.}
\item{num.cores}{an optional numeric. The number of processes executed in parallel if 
\code{parallel} is \code{TRUE}. By default it is set equal to the number of CPU cores.}
}

\details{
This function is based on the discussion given in López-de-Lacalle (2013)
about the Expectation-Maximization (EM) algorithm in the context of structural time series models
The previous reference includes an \R package called \pkg{stsm.em} but since it 
uses some non-standard options it has not been distributed.
A version of the most relevant procedures is provided in the function \code{maxlik.em}.

The traditional design of the EM algorithm in the context of structural 
time series models is run with \code{type = "standard"}.
This approach is described for instance in Durbin and Koopman (2001) Section 7.3.4.
A modified version introduced in López-de-Lacalle (2013) can be run 
using \code{type = "modified"}. A mixture of both approaches is also possible by 
setting \code{type = "mix"}. In that case, the modified version is run at those 
iterations indicated in \code{mod.steps} and the traditional version is run 
in the remaining iterations. As we do not know beforehand the number of iterations 
required for convergence, \code{mod.steps} should be defined considering up to 
\code{maxiter} possible iterations.

In pure variance models, the Kalman filter may converge to a steady state.
The parameters in \code{kfconv} control how convergence of the filter is determined. 
It is considered that convergence is reached when the following is observed:
the change in the variance of the prediction error over the last \code{kfconv[2]} 
consecutive iterations of the filter is below the tolerance \code{kfconv[1]}.
The iteration at which the the Kalman smoother has converged is the 
iteration where the Kalman filter converged multiplied by the factor 
\code{kfconv[3]}. If provided, \code{kfconv[3]} should be equal or greater than unity.

The argument \code{r.interval} is a three-length vector.
The first two elements are the initial lower and upper ends of the interval where the 
variance parameters are searched. 
The third element is the first iteration of the EM algorithm after which the initial 
searching interval is narrowed. A relatively width initial interval is recommended.
As the algorithm makes progress, the interval is automatically narrowed according to 
the values and path followed in the first \code{r.interval[3]} and subsequent iterations.
}

\value{
A list containing the elements:
\item{Mpars}{a matrix with the values of the parameters stored by row for each iteration 
of the procedure.}
\item{pars}{parameter values at the local optimum or point where the algorithm stopped.}
\item{iter}{number of iterations until convergence or stopping.}
%
For \code{type = "modified"} or \code{type = "mix"} the element
\code{calls.v1} is also returned. It reports, for each parameter, the number of calls
that were done to the traditional version due to failure to convergence of the 
root finding procedure.
}

\seealso{
\code{\link{maxlik.td}}.
}

\references{
Durbin, J. and Koopman, S.J. (2001). Section 7.3.4.
\emph{Time Series Analysis by State Space Methods}. 
Oxford University Press.

Harvey, A. C. (1989). Section 4.2.4.
\emph{Forecasting, Structural Time Series Models and the Kalman Filter}. 
Cambridge University Press.

Koopman, S.J. and Shephard, N. (1992)
Exact Score for Time Series Models in State Space Form.
\emph{Biometrika} \bold{79}(4), pp. 823-826.

Koopman, S.J. (1993). 
Disturbance Smoother for State Space Models.
\emph{Biometrika} \bold{80}(1), pp. 117-126.

López-de-Lacalle, J. (2013).
\sQuote{Why Does the Expectation-Maximization Algorithm Converge Slowly 
in Pure Variance Structural Time Series Models?} Unpublished manuscript.

Shumway, R.H. and Stoffer, D.S. (1982)
An Approach to Time Series Smoothing and Forecasting Using the EM Algorithm.
\emph{Journal of Time Series Analysis} \bold{3}(4), pp. 253-264.
}

\examples{
# fit a local level plus seasonal model to a simulated sample series
# using the three versions of the EM algorithm
# the same solution is found by all versions (up to a tolerance)
# the modified version converges in fewer iterations, yet it involves 
# more computations
data("llmseas")
m <- stsm.class::stsm.model(model = "llm+seas", y = llmseas, 
  pars = c(var1 = 1, var2 = 1, var3 = 1))

# original version
res1 <- maxlik.em(m, type = "standard", 
  tol = 0.001, maxiter = 350, kfconv = c(0.01, 10, 1))
res1$pars
res1$iter

# modified version
res2 <- maxlik.em(m, type = "modified", 
  tol = 0.001, maxiter = 250, kfconv = c(0.01, 10, 1),  
  ur.maxiter = 1000, r.interval = c(0.001, var(m@y)*0.75, 20))
res2$pars
res2$iter
res2$calls.v1

# mixture, the modified version is run every 10 iterations
# starting in the third one
res3 <- maxlik.em(m, type = "mix", 
  tol = 0.001, maxiter = 250, kfconv = c(0.01, 10, 1),  
  ur.maxiter = 1000, r.interval = c(0.001, var(m@y)*0.75, 20),
  mod.steps = seq(3, 200, 10))
res3$pars
res3$iter
res3$calls.v1
}

\keyword{ts}
\keyword{models}
\keyword{nonlinear}
