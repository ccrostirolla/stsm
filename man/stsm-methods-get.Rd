\name{stsm-get-methods}
\docType{methods}
\alias{stsm-get-methods}
\alias{get.cpar}
\alias{get.cpar,stsm-method}
\alias{get.nopars}
\alias{get.nopars,stsm-method}
\alias{get.pars}
\alias{get.pars,stsm-method}

\title{Getter Methods for Class \code{stsm}}

\description{
Get access to the information stored in the slots \code{cpar}, \code{nopars} and 
\code{pars} in objects of class \code{\link{stsm}}.
}

\usage{
\S4method{get.cpar}{stsm}(x, rescale = FALSE)
\S4method{get.nopars}{stsm}(x, rescale = FALSE)
\S4method{get.pars}{stsm}(x, rescale = FALSE, gradient = FALSE)
}

\arguments{
\item{x}{an object of class \code{\link{stsm}}.}
\item{rescale}{logical. If \code{TRUE}, relative variance parameters are 
rescaled into absolute variances. Ignored if \code{x@cpar} is null.}
\item{gradient}{logical. If \code{TRUE}, first order derivatives of 
\code{\link{transPars}} with respect to the parameters in the slot \code{pars} are returned.}
}

\details{
\emph{Transformation of the parameters of the model.}
The method \code{\link{transPars}} allows parameterizing the model 
in terms of an auxiliary vector of parameters.
The output of \code{get.pars} is returned in terms of the actual parameters of the model, i.e.,
the variances and the autoregressive coefficients if they are part of the model.
With the standard parameterization, \code{x@transPars = NULL}, \code{get.pars(x)} returns 
the output stored in \code{x@pars}.
When the model is parameterized in terms of an auxiliary set of parameters 
\eqn{\theta}, \code{get.pars} return the variance parameters instead of the values of 
\eqn{\theta} that are stored in \code{x@pars}.
For example, with \code{x@transPars = "square"} (where the variances are \eqn{\theta^2}),
\code{ger.pars} returns \eqn{\theta^2} while \code{x@pars} contains the vector \eqn{\theta}.

\emph{Absolute and relative variances.}
%NOTE this paragraph is included in 'sts-model.Rd', 
%repeated here for introduction and a standalone explanation
The model can be defined in terms of relative variances.
In this case, the variance that acts as a scaling parameter is stored in the slot \code{cpar}. 
Otherwise, \code{cpar} is null and ignored.
Typically, the scaling parameter will be chosen to be the variance parameter that is 
concentrated out of the likelihood function.
%
If \code{rescale = TRUE}, the relative variance parameters are rescaled into 
absolute variance parameters (i.e., they are multiplied by \code{x@cpar}) 
and then returned by these methods. 
If \code{rescale = FALSE}, relative variance parameters are returned, that is,
the variances divided by the scaling parameter \code{cpar}.
Since the scaling parameter is one of the variances, the relative variance stored 
in \code{cpar} is \eqn{1} (the parameter divided by itself).

\emph{Transformation of parameters in a model defined in terms of relative 
variances.} When a model is defined so that the parameters are the relative variances 
(\code{cpar} is not null) and a parameterization \code{transPars} is also specified,
then the transformation of parameters is applied to the relative variances,
not to the absolute variances. The relative variances  are first transformed and 
afterwards they are rescaled back to absolute variances if requested by 
setting \code{rescale = TRUE}. 
The transformation \code{transPars} is applied to the parameters defined in 
\code{pars}; \code{cpar} is assumed to be chosen following other rationale; usually,
it is the value that maximizes the likelihood since one of the variance parameters 
can be concentrated out of the likelihood function.

\emph{Note.}
When \code{cpar} is not null, it is more convenient to store in the slots 
\code{pars} and \code{nopars} the values of the relative variances, while
the slot \code{cpar} stores the value of the scaling parameter rather than 
the relative variance (which will be \eqn{1}).
If the relative values were stored, then the scaling parameter would 
need to be recomputed each time the value is requested by \code{get.cpar}.
Assuming that \code{cpar} is the parameter that is concentrated out of the likelihood function,
the expression that maximizes the likelihood should be evaluated whenever 
the value is requested to be printed or to do any other operation.
To avoid this, the scaling value is directly stored. This approach makes also 
sense with the way the method \code{\link{set.cpar}} works.

\emph{Note for users.}
For those users that are not familiar with the design and internal structure of the 
class \code{\link{stsm}}, it is safer to use the get and set methods 
rather than retrieving or modifying the contents of the slots through the \code{@} and 
\code{@<-} operators.
}

\value{
\item{get.cpar}{named numeric of length one.}
\item{get.nopars}{named numeric vector.}
\item{get.pars}{named numeric vector.}
}

\seealso{
\code{\link{stsm}}.
}

\examples{
# sample models with arbitrary parameter values

# model in standard parameterization
# internal parameter values are the same as the model parameter
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 15, "var3" = 30))
m@pars
get.pars(m)

# model parameterized, the variances are the square
# of an auxiliary vector of parameters
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 15, "var3" = 30), transPars = "square")
# auxiliary vector of parameters
m@pars
# parameters of the model, variances
get.pars(m)

# model rescaled, variances are relative to 'var1'
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var2" = 15, "var3" = 30), cpar = c("var1" = 2))
# internal values
m@pars
m@cpar
# relative variances
get.pars(m)
get.cpar(m)
# absolute variances
get.pars(m, rescale = TRUE)
get.cpar(m, rescale = TRUE)

# model defined in terms of relative variances 
# and with the parameterization \code{transPars="square};
# the transformation is applied to the relative variances,
# the relative variances are first transformed and afterwards 
# they are rescaled back to absolute variances if requested
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var2" = 3, "var3" = 6), cpar = c("var1" = 2),
  transPars = "square")
c(get.cpar(m, rescale = FALSE), get.pars(m, rescale = FALSE))
c(get.cpar(m, rescale = TRUE), get.pars(m, rescale = TRUE))

# when 'cpar' is defined, 'nopars' is also interpreted as a relative variance
# and therefore it is rescaled if absolute variances are requested
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var2" = 3), cpar = c("var1" = 2), nopars = c("var3" = 6),
  transPars = NULL)
v <- c(get.cpar(m, rescale = FALSE), get.pars(m, rescale = FALSE), get.nopars(m, rescale = FALSE))
v[c("var1", "var2", "var3")]
v <- c(get.cpar(m, rescale = TRUE), get.pars(m, rescale = TRUE), get.nopars(m, rescale = TRUE))
v[c("var1", "var2", "var3")]

# 'nopars' is rescaled as shown in the previous example
# but it is not affected by the parameterization chosen for 'pars'
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var2" = 3), cpar = c("var1" = 2), nopars = c("var3" = 6),
  transPars = "square")
v <- c(get.cpar(m, rescale = FALSE), get.pars(m, rescale = FALSE), get.nopars(m, rescale = FALSE))
v[c("var1", "var2", "var3")]
v <- c(get.cpar(m, rescale = TRUE), get.pars(m, rescale = TRUE), get.nopars(m, rescale = TRUE))
v[c("var1", "var2", "var3")]
}

\keyword{methods}
