\name{stsm-transPars-methods}
\docType{methods}
\alias{stsm-transPars-methods}
\alias{transPars,generic-method}
\alias{transPars,numeric-method}
\alias{transPars,stsm-method}
\alias{transPars}

\title{Parameterization of Models Defined in the Class \code{stsm}}

\description{
This method provides different transformations of the parameters of 
a structural time series model. 
}

\usage{
\S4method{transPars}{generic}(x, 
    type = c("square", "StructTS", "exp", "exp10sq"),
    gradient = FALSE, hessian = FALSE, 
    rp, sclrho = 1.7, sclomega = 1.7, ftrans = NULL, ...)
\S4method{transPars}{numeric}(x, 
    type = eval(formals(transPars)$type),
    gradient = FALSE, hessian = FALSE, 
    rp, sclrho = 1.7, sclomega = 1.7, ftrans = NULL, ...)
\S4method{transPars}{stsm}(x, type = NULL,
    gradient = FALSE, hessian = FALSE, 
    rp, sclrho = 1.7, sclomega = 1.7, ftrans = NULL, ...)
}

\arguments{
\item{x}{an object of class \code{\link{stsm}}.}
\item{type}{a character string indicating the type of transformation. 
Ignored if \code{x} is of class \code{stsm}. See details.}
\item{gradient}{logical. If \code{TRUE}, first order derivatives of the 
transformation function with respect to the parameters in the slot \code{pars} are returned.}
\item{hessian}{logical. If \code{TRUE}, second order derivatives of the 
transformation function with respect to the parameters in the slot \code{pars} is returned.}
\item{rp}{numeric value. Regularization parameter used with \code{type = StrucTS}. 
By default it is the variance of the data \code{x@y} divided by \eqn{100}.}
\item{sclrho}{numeric value. Currently ignored.}
\item{sclomega}{numeric value. Currently ignored.}
\item{ftrans}{a function defining an alternative transformation of the parameters. 
Ignored if \code{x} is of class \code{stsm}.
See example below.}
\item{...}{additional arguments to be passed to \code{ftrans}.}
}

\details{
Rather than using the standard parameterization of the model
(in terms of variances and autoregressive coefficients if they are part of the model),
it can be parameterized in terms of an auxiliary set of parameters \eqn{\theta}.
This may be useful for example when the parameters of the model 
are selected by means of a numerical optimization algorithm. Choosing a suitable 
parameterization ensures that the solution returned by the algorithm meets
some constraints such as positive variances or autoregressive coefficients 
within the region of stationarity.

The method \code{transPars} can be applied both on a named vector of parameters,
e.g. \code{x@pars} or on a model of class \code{\link{stsm}}.

When the slot \code{transPars} is not null, the model is parameterized in 
terms of \eqn{\theta}. The following transformation of 
parameters can be considered:
\itemize{
\item{\code{"square"}: the variance parameters are the square of \eqn{\theta}.}
\item{\code{"StructTS"}: transformation used in the function \code{\link[stats]{StructTS}} 
of the \pkg{stats} package.}
\item{\code{"exp"}: the variance parameters are the exponential of \eqn{\theta}.}
\item{\code{"exp10sq"}: the variance parameters are \eqn{(exp(-\theta)/10)^2}.}
}

In the model \code{trend+ar2} defined in \code{\link{stsm.model}},
the autoregressive coefficients, \eqn{\phi}, are transformed to lie in the region of 
stationarity: given \eqn{z1 = \phi_1 / (1 + |\phi_1|)}, \eqn{z2 = \phi_2 / (1 + |\phi_2|)},
the transformed coefficients are \eqn{\phi_1^* = z1 + z2} and 
\eqn{\phi_2 = - z1 \cdot z2}.

Other transformations can be defined through the argument \code{ftrans}, which can also 
be defined in the slot \code{transPars} of a \code{stsm} object.
\code{ftrans} must be a function returning a list containing an element 
called \code{pars} and two other optional elements called \code{gradient} and 
\code{hessian}. The parameters to be transformed are identified by their names.
The variances follow the naming convention of the regular expression
\dQuote{^var\\d\{1,2\}$}, e.g. \code{var1}, \code{var2},... 
The variances of the initial state vector may also be transformed 
if they are included in the slot \code{pars}; their names follow a similar 
naming convention, \code{P01}, \code{P02},...
An example of \code{ftrans} is given below.

\emph{Note:} If a transformation is defined by means of \code{ftrans}
the user may need to update the slots \code{lower} and \code{upper}
if some bounds are still applied to the auxiliary parameters. For example,
\code{transPars="StructTS"} does not always yield positive variances and
hence lower bounds equal to \code{0} are needed. By default lower and upper
bounds are not considered if \code{ftrans} is used.
%(they are set equal to minus infinity and infinity)

The output of \code{\link{get.pars}} is given in terms of 
the actual parameters of the model.
For example, if the model is parameterized so that \eqn{\theta^2} are the
variances of the model and \eqn{\theta} are the auxiliary parameters then, 
the slot \code{pars} contains the values of \eqn{\theta} and 
\code{ger.pars} returns \eqn{\theta^2}.

The transformation \code{transPars} is applied to the parameters included
in the slot \code{pars}. The transformation does not affect \code{nopars} and \code{cpar}.
The former slot is considered fixed while the latter will in practice be set equal 
to a particular value, for example the value that maximizes the concentrated 
likelihood function, for which a specific expression can be obtained.
}

\value{
A list containing a named numeric vector with the values of the transformed parameters.
If requested, the gradient and Hessian of the transformation function with 
respect to the parameters are returned.
}

\seealso{
\code{\link{stsm}},
\code{\link{get.pars}}.
}

\examples{
# sample models with arbitrary parameter values

# model in standard parameterization
# ower bounds imposed on the variance parameters
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 15, "var3" = 30), transPars = NULL)
get.pars(m)
m@lower

# square transformation
# negative values are allowed in 'pars' since 
# the square will yield positive variances
# in fact no lower bounds need to be imposed on the auxiliary parameters
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = -2, "var2" = -5, "var3" = 10), transPars = "square")
validObject(m)
m@lower
m@pars
get.pars(m)

# 'ftrans', alternative transformation of parameters;
# the following parameterization is sometimes found:
# variance = exp(-theta) / 10
# the function 'ftrans' following the rules given in the details 
# above can be defined as follows:

ftrans <- function(x, gradient = FALSE, hessian = FALSE)
{
  tpars <- x
  p <- length(x)
  nmspars <- names(x)
  idvar <- grep("^var|P0\\\\d{1,2}$", nmspars, value = FALSE)

  if (gradient) {
    d1 <- rep(NA, p)
    names(d1) <- nmspars
  } else d1 <- NULL
  if (hessian) {
    d2 <- matrix(0, p, p)
    rownames(d2) <- colnames(d2) <- nmspars
  } else d2 <- NULL

  if (length(idvar) > 0) {
    tpars[idvar] <- exp(-x[idvar]) / 10
  } else warning("No changes done by 'transPars'.")

  if (gradient)
  {
    if (length(idvar) > 0)
      d1[idvar] <- -tpars[idvar]
  }
  if (hessian) {
    diag(d2)[idvar] <- tpars[idvar]
  }
  
  list(pars = tpars, gradient = d1, hessian = d2)
}

# now 'ftrans' can be passed to 'transPars' and be applied
# on a named vector of parameters or on a 'stsm' object
transPars(c("var1" = 2, "var2" = 15, "var3" = 30), 
  ftrans = ftrans, gradient = TRUE, hessian = TRUE)
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 15, "var3" = 30), transPars = ftrans)
get.pars(m)
}

\keyword{methods}
