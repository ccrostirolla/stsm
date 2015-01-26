\name{stsm-validObject-methods}
\docType{methods}
\alias{stsm-validObject-methods}
\alias{check.bounds}
\alias{check.bounds,stsm-method}
\alias{validObject}
\alias{validObject,stsm-method}

\title{Check the Validity of an Object of Class \code{stsm}}

\description{
Methods to check the validity of an object of class \code{\link{stsm}}.
}

\usage{
\S4method{check.bounds}{stsm}(x)
\S4method{validObject}{stsm}(object)
}

\arguments{
\item{x}{an object of class \code{\link{stsm}}.}
\item{object}{an object of class \code{\link{stsm}}.}
}

\details{
\code{check.bounds} checks that the values in the slot \code{pars}
lie within the lower and upper bounds. These bounds are stored in the slots 
\code{lower} and \code{upper}. Default values or specific values can be given 
when creating the object by means of \code{\link{stsm.model}}.

\code{check.bounds} is called by \code{validObject}. 
In some settings it may be required to check only that the parameters 
are within the required bounds.

\code{validObject} checks additional requirements: e.g. all the parameters 
taking part in the selected model are either in the slots \code{pars}, 
\code{nopars} or \code{cpar}; 
%(parameters are identified by their names in the slots)
it is also checked that the parameters are no duplicated in those slots.

This method is called by \code{\link{stsm-set-methods}} defined 
for the slots \code{pars}, \code{nopars} or \code{cpar}. That's why it is safer 
to use the setter methods instead of a direct modification through the operator \code{@<-}.
}

\value{
If the input object is valid according to the class definition, 
the logical \code{TRUE} is returned. Otherwise, an error message is 
returned.
}

\seealso{
\code{\link{stsm}} and
examples in \code{\link{stsm-set-methods}}.
}

\examples{
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 15, "var3" = 30))
validObject(m)

# force a wrong value (negative variance)
m@pars[1] <- -1
try(validObject(m))
try(check.bounds(m))

# duplicates not allowed
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 15, "var3" = 30))
# try to define 'var1', already in 'pars', in the slot 'nopars'
try(m <- set.nopars(m, c(var1=22)))
# force a duplicate
m@nopars <- c(m@nopars, var1 = 22)
try(validObject(m))
}

\keyword{methods}
