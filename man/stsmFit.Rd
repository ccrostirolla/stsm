\name{stsmFit}
\alias{stsmFit}

\title{Interface to Different Fitting Procedures}

\description{
This function is an interface to the methods available in the 
package for fitting a structural time series model (maximum likelihood 
in the time and frequency domain via different optimization algorithms). 
}

\usage{
stsmFit(x, stsm.method = c("maxlik.fd.scoring", "maxlik.td.scoring", 
  "maxlik.fd.optim", "maxlik.td.optim"), xreg = NULL, ...)
}

\arguments{
\item{x}{an object of class \link{stsm}.}
\item{stsm.method}{a character indicating the method to be used.}
\item{...}{further arguments to be passed to the the function selected in 
\code{stsm.method}.}
\item{xreg}{vector or matrix of external regresors.}
}

\details{
This interface is useful to simplify the code and reduce the number of 
arguments in functions that call some of those functions that can be specified 
through \code{stsm.method}. 
For example, the package \pkg{tsoutliers} uses:
\code{do.call("stsmFit", args = c(list(x = y, args.tsmethod))}, 
where \code{args.tsmethod} is a list containing the arguments to 
be passed to \code{stsmFit}, which includes \code{stsm.method}.
Thus, the code is simplified since no \code{switch} or \code{if} 
statements are necessary; the number of arguments is also reduced 
since those that are passed tot the function specified in 
\code{stsm.method} are gathered in a list.

The external regressors can be defined in the input object 
\code{x} of class \link{stsm}. This is the way the recommend specification 
for functions \link{maxlik.td} and \link{maxlik.fd}. This interface allows 
defining \code{xreg} as an argument passed to this function because it 
simplifies the code in some functions of package \pkg{tsoutliers}. 
If \code{xreg} and \code{x@xreg} are both not \code{NULL} and error is returned.
}

\value{
A list of class \code{stsmFit}. 
See the section \sQuote{Value} in \code{\link{maxlik.fd}}.
}

\seealso{
\code{\link{maxlik.fd.scoring}},
\code{\link{maxlik.td.scoring}},
\code{\link{maxlik.fd.optim}} and
\code{\link{maxlik.td.optim}}.
}

\keyword{ts}
\keyword{models}
