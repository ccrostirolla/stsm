\name{methods-stsmFit}
\alias{methods-stsmFit}
\alias{coef.stsmFit}
\alias{print.stsmFit}
\alias{fitted.stsm}
\alias{fitted.stsmFit}
\alias{residuals.stsmFit}
\alias{plot.stsmComponents}
\alias{predict.stsm}
\alias{predict.stsmFit}
\alias{plot.stsmPredict}
\alias{tsSmooth.stsm}
\alias{tsSmooth.stsmFit}
\alias{plot.stsmSmooth}
\alias{tsdiag.stsmFit}

\title{Methods to Extract Information from a Fitted \code{stsm} Model Object}

\description{
Common methods to print and display results for an object of class 
\code{\link{stsm}} or \code{stsmFit} returned by 
functions \code{\link{maxlik.fd}} and \code{\link{maxlik.fd}}.
}

\usage{
\method{coef}{stsmFit}(object, ...)
\method{print}{stsmFit}(x, digits = max(3L, getOption("digits") - 3L), 
  vcov.type = c("hessian", "infomat", "OPG", "sandwich", "optimHessian"), ...)
\method{fitted}{stsmFit}(object, std.rediduals = TRUE, version = c("KFKSDS", "stats"), ...)
\method{residuals}{stsmFit}(object, standardised = FALSE, version = c("KFKSDS", "stats"), ...)
\method{plot}{stsmComponents}(x, ...)
\method{predict}{stsm}(object, n.ahead = 1L, se.fit = TRUE, 
  version = c("KFKSDS", "stats"), ...)
\method{predict}{stsmFit}(object, n.ahead = 1L, se.fit = TRUE, 
  version = c("KFKSDS", "stats"), ...)
\method{plot}{stsmPredict}(x, ...)
\method{tsSmooth}{stsm}(object, version = c("KFKSDS", "stats"), ...)
\method{tsSmooth}{stsmFit}(object, version = c("KFKSDS", "stats"), ...)
\method{plot}{stsmSmooth}(x, ...)
\method{tsdiag}{stsmFit}(object, gof.lag = 10L, ...)
}

\arguments{
\item{object}{an object of class \code{\link{stsm}} or a \code{stsmFit} list.}
\item{x}{a \code{stsmFit} list. For \code{plot} methods it is a list of class 
\code{stsmComponents}, \code{stsmPredict} or \code{stsmSmooth}.}
\item{digits}{minimal number of significant digits, see \code{print.default}.}
\item{vcov.type}{a character indicating the type of covariance matrix to be used to 
  compute the standard errors of the parameter estimates.}
\item{version}{a character indicating whether the Kalman filter and smoother functions 
from package \pkg{KFKSDS} or from the \pkg{stats} package should be used.}
\item{std.rediduals}{logical. If \code{TRUE} residuals are standardised.}
\item{standardised}{logical. If \code{TRUE} standardised are returned.}
\item{n.ahead}{a numeric, number of observations ahead to perform prediction.}
\item{se.fit}{logical. If \code{TRUE} standard errors of the predictions are returned.}
\item{gof.lag}{numeric, number of lag autocorrelation coefficients to which apply the Box test.}
\item{type}{A character. Type of information used to compute the covariance matrix of the 
parameters of the fitted model: information matrix, Hessian or a mixture of them.}
\item{...}{additional arguments to be passed to the functions called in these methods.}
}

\details{
These methods are based on those with the same name available for the output
returned by \code{\link[stats]{StructTS}} in the \pkg{stats} package.

%Some of these methods provide information upon the \code{stsmFit} list returned by 
%\code{\link{maxlik.fd}} and \code{\link{maxlik.td}}. 
%The data and results are stored in the 
%\code{\link{stsm}} object of the output list, which is the 
%model updated with the optimal parameter values. Hence, these methods 
%can be either applied on a \code{\link{stsm}} object or 
%a \code{stsmFit} list.
%
These methods are originally intended to provide summary information from 
a model fitted by maximum likelihood. 
Thus, the most natural input for them is a \code{stsmFit} list returned by 
\code{\link{maxlik.fd}} and \code{\link{maxlik.td}}.
Nevertheless, as the information and the data required by these methods 
are available in the slots of a \code{\link{stsm}} object, they 
can also be applied directly on an object of class \code{stsm}. 
This can be useful, for example, when we know a set of parameter estimates 
that was obtained from a method other than \code{maxlik.fd} or \code{maxlik.td}. 
By simply updating the slot \code{pars} of the \code{stsm} object, 
the residuals and the filtered and smoothed components
are readily available through these methods.

In most cases, the ellipsis, \code{...}, is kept in the definitions of the methods 
just because it is part of the parent method. It has some functionality nonetheless.
For the methods \code{fitted.stsm} and \code{predict.stsm} it can be used to pass 
argument \code{P0cov} to function \code{\link{char2numeric}}. 
It can also be used to pass graphical parameters to \code{par} in 
method \code{plot.stsmComponents} and \code{plot.stsmSmooth} or to \code{plot}
in \code{plot.stsmPredict} and \code{plot.stsmSmooth}.
%For \code{coef.stsmFit} and \code{print.stsmFit} it can be used to pass 
%argument \code{rescale} to \code{\link{get.pars}}.

By default in method \code{fitted}, \code{std.rediduals = TRUE} so that 
it behaves as in previous versions of the package where this 
argument was not available.
In method \code{rediduals}, by default \code{standardised = FALSE} 
because it is more convenient when used in package \pkg{tsoutliers}
(this argument does not need to be explicitly defined and hence the 
method residuals is used in the same way both on \code{\link[stats]{arima}} 
models and for \code{stsm}).
}

\value{
The following information is returned by these methods:
%
\item{coef.stsmFit, print.stsmFit}{optimal parameter values.}
\item{fitted.stsm, fitted.stsmFit}{an list of class \code{stsmComponents} 
containing the filtered components and residuals.}
\item{residuals.stsmFit}{residuals in the fitted model.}
\item{plot.stsmComponents}{plot of the filtered components.}
\item{predict.stsm, predict.stsmFit}{predictions of the input time series and standard errors.}
\item{plot.stsmPredict}{plot of the predictions.}
\item{tsSmooth.stsm, tsSmooth.stsmFit}{smoothed components.}
\item{plot.stsmSmooth}{plot of the smoothed components.}
\item{tsdiag.stsm, tsdiag.stsmFit}{plot of diagnostic tests.}
}

\seealso{
\code{\link{confint.stsmFit}},
\code{\link{maxlik.fd}},
\code{\link{maxlik.td}},
\code{\link{vcov.stsmFit}},
\code{\link{stsm}}.
}

\examples{
# fit the local level plus seasonal model to a 
# sample simulated series
data("llmseas")
m <- stsm.model(model = "llm+seas", y = llmseas)
res <- maxlik.fd.scoring(m = m, step = NULL, 
  information = "expected", control = list(maxit = 100, tol = 0.001))
print(res)

#diagnostic
tsdiag(res)

# display estimated components with 95% confidence bands
comps <- tsSmooth(res)
plot(comps)
title(main = "smoothed trend and seasonal components")

# plot predictions eight periods ahead
pred <- predict(res, n.ahead = 8, se.fit = TRUE)
plot(pred)
}

\keyword{ts}
