\name{stsm-show-methods}
\docType{methods}
\alias{stsm-show-methods}
\alias{show}
\alias{show,stsm-method}

\title{Display an Object of Class \code{stsm}}

\description{
This method displays summary information about 
an object of class \code{\link{stsm}}.
}

\usage{
\S4method{show}{stsm}(object)
}

\arguments{
\item{object}{an object of class \code{\link{stsm}}.}
}

\details{
A succinct summary of the object (name of the model and parameter values)
is printed.
}

\value{
Invisible \code{NULL}.
}

\seealso{
\code{\link{stsm-class}}.
}

\examples{
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 15, "var3" = 30))
show(m)
# or just
m
}

\keyword{methods}

