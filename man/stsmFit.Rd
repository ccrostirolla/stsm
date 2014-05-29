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
  "maxlik.fd.optim", "maxlik.td.optim"), ...)
}

\arguments{
\item{x}{a time series.}
\item{stsm.method}{a character indicating the method to be used.}
\item{...}{further arguments to be passed to the the function selected in 
\code{stsm.method}.}
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
