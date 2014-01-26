\name{methods-vcov-confint}
\alias{methods-vcov-confint}
\alias{vcov.stsmFit}
\alias{vcov.stsm}
\alias{confint.stsmFit}

\title{Variance-covariance Matrix for a Fitted \code{stsm} Model Object}

\description{
The method \code{vcov} computes the variance-covariance matrix of the parameters fitted 
in a structural time series model. This matrix is used to compute confidence 
intervals for those parameters returned by the \code{\link[stsm]{coef.stsmFit}} method.
}

\usage{
\method{vcov}{stsmFit}(object, 
  type = c("hessian", "infomat", "OPG", "sandwich", "optimHessian"), ...)
\method{vcov}{stsm}(object, 
  type = c("hessian", "infomat", "OPG", "sandwich"), 
  domain = c("frequency", "time"), ...)
\method{confint}{stsmFit}(object, parm, level = 0.95, 
  type = c("vcov", "bootstrap"),
  vcov.type = c("hessian", "infomat", "OPG", "sandwich", "optimHessian"), 
  breps = 100, ...)
}

\arguments{
\item{object}{a \code{stsmFit} list or object of class \code{\link[stsm.class]{stsm}}.}
%
\item{type}{a character. In \code{vcov}, it is the type of covariance matrix.
In \code{confit}, it is the type of confidence intervals: based on the 
covariance matrix of the estimated parameters or on a bootstrap procedure.}
%
\item{domain}{a character indicating whether the covariance matrix is obtained upon the 
frequency or time domain likelihood function.}
%
\item{parm}{character indicating the name of the parameter on to obtain the confidence interval.
If missing, all parameters are considered.}
%
\item{level}{the confidence level.}
%
\item{vcov.type}{a character indicating the type of covariance matrix.
Ignored if \code{type = "bootstrap"}.}
%
\item{breps}{number of bootstrap replicates. Ignored if \code{type "vcov"}.}
%
\item{...}{additional arguments to be passed to the functions called in these methods.
Currently ignored.}
}

\details{
The following estimators of the covariance matrix of parameter estimates 
are available (Davidson and MacKinnon (2004), Section 10.4):
%
\itemize{
\item{\code{hessian}: the inverse of the analytical Hessian.}
%
\item{\code{infomat}: the inverse of the analytical expression for the information matrix.}
%
\item{\code{OPG}: the inverse of the outer product of the analytical gradient.
Also known as the BHHH estimator since it was proposed by 
Berndt, Hall, Hall and Hausman (1974).
This method requires only first order derivatives. It tends to be less
reliable in small samples.}
%
\item{\code{sandwich}: the sandwich estimator defined as:
\eqn{H^{-1}  (G'G)  H^{-1},}
where \eqn{G} is the gradient vector and \eqn{H} is the Hessian matrix.
It requires more computations and
may be unreliable in small samples. However, contrary to the previous methods,
it is valid when the information matrix equality does not hold
for example due to misspecification of the model.}
%
\item{\code{optimHessian}: the inverse of the numerical Hessian 
returned by \code{\link[stats]{optim}}.}
}

The natural input for the method method \code{vcov} is a \code{stsmFit} 
list returned by \code{\link{maxlik.fd}} or \code{\link{maxlik.td}}.
However, \code{vcov} can be also applied directly on a 
\code{\link[stsm.class]{stsm}} model object.
This is useful for example when maximum likelihood parameter estimates 
are found by means of \code{\link[stats]{optim}} using the functions 
\code{maxlik.fd.optim} or \code{maxlik.fd.optim}. In that case,
only the covariance matrix based on the numerical Hessian returned by 
\code{optim} would be available. Updating the slot \code{pars} of a 
\code{stsm} model object with the parameter values obtained from other 
algorithm is a convenient solution to obtain the covariance matrix 
based on other methods available in \code{vcov.stsmFit}.

For the time domain likelihood function the covariance matrix of the 
initial state vector is considered diagonal, \code{P0cov = FALSE}.

The analytical Hessian for the time domain version is not available, 
the information matrix is used instead.

By default, \code{vcov.type = "infomat"} for the time domain likelihood function 
and \code{vcov.type = "hessian"} for the frequency domain likelihood function.

Confidence intervals can either be computed upon the covariance 
matrix of the parameter estimates (\code{type = "vcov"}) or by means of 
bootstrapping (\code{type = "bootstrap"}).
The bootstrap approach takes advantage of the following result 
(Harvey (1989) eq. (4.3.25)):
%
\deqn{4\pi I(\lambda_j) / g(\lambda_j) \sim \chi^2_2, 
  \hbox{ for } j\neq 0, n/2 \hbox{ (for n even)}}
\deqn{2\pi I(\lambda_j) / g(\lambda_j) \sim \chi^2_1, 
  \hbox{ for } j=0, n/2 \hbox{ (for n even)}}
%
where \eqn{I(\lambda_j)} and \eqn{g(\lambda_j)} are respectively the periodogram and 
the spectral generating function at frequency \eqn{\lambda_j}.
Upon this result, bootstrap replicates of the periodogram are generated and 
for each of them parameter estimates are obtained maximizing the spectral 
likelihood function. The quantiles of the bootstrapped parameter estimates are 
the confidence interval.
Dahlhaus and Janas (1996) studied the properties of the frequency domain bootstrap
which has been applied, among others, in Koopman and Wong (2006).
An advantage of the bootstrap method is that it yields confidence intervals
within the bounds of the parameters, i.e., positive variances.
This procedure is computationally intensive and requires some time to run,
especially for large \code{breps}.

%A barrier term is not considered since it is an element added to the likelihood
%to constrain the searching parameter space, given the local optimum found inside 
%the feasible area the barrier is no longer needed to evaluate the likelihood,
%thus the barrier is not considered here to avoid that if distort the evaluation of 
%the covariance matrix and the standard errors.
}

\value{
\item{vcov.stsm, vcov.stsmFit}{return the covariance matrix of the parameters of the model.}
\item{confint.stsmFit}{returns a matrix containing confidence intervals 
for the parameters of the model.}
}

\seealso{
\code{\link{maxlik.fd}},
\code{\link{maxlik.td}},
\code{\link{methods-stsmFit}},
\code{\link[stsm.class]{stsm-class}}.
}

\references{
Berndt, E. R., Hall, B. H., Hall, R. E. and Hausman, J. A. (1974).
\sQuote{Estimation and inference in nonlinear structural models}.
\emph{Annals of Economic and Social Measurement}, \bold{3}, pp. 653-65.

Dahlhaus, R. and Janas, D. (1996).
\sQuote{A Frequency Domain Bootstrap for Ratio Statistics in Time Series Analysis}.
\emph{Annals of Statistics}, \bold{24}(5), pp. 1934-1963.

Davidson, R. and MacKinnon, J. G. (2004). Section 10.4.
\emph{Econometric Theory and Methods}. Oxford University Press.

Koopman, S. J. and Wong, S. Y. (2006).
\sQuote{Extracting Business Cycles using Semi-Parametric Time-varying Spectra with 
Applications to US Macroeconomic Time Series}. 
Tinbergen Institute Discussion Papers, No. 2006-105/4.
\url{http://papers.tinbergen.nl/06105.pdf}
}

\examples{
\dontrun{
data("llmseas")
# fit the local level plus seasonal model with default arguments
# using the Newton-Raphson algorithm
m <- stsm.model(model = "llm+seas", y = llmseas[,1])
res <- maxlik.fd.scoring(m = m, information = "observed")
coef(res)
# confidence intervals for parameter estimates ...
# ... based on the covariance matrix of parameter estimates
# gives a warning since the lower limit of the confidence interval
# for parameter 'var2' was forced to be non-negative (fixed to 0)
civcov <- confint(res, type = "vcov", vcov.type = "hessian")
civcov
# ... based on bootstrapping the periodogram
# NOTE: this will take a while to run
set.seed(643)
ciboot <- confint(res, type = "bootstrap", breps = 100)
ciboot
}
}

\keyword{ts}
