\name{stsm-set-methods}
\docType{methods}
\alias{stsm-set-methods}
\alias{set.cpar}
\alias{set.cpar,stsm-method}
\alias{set.nopars}
\alias{set.nopars,stsm-method}
\alias{set.pars}
\alias{set.pars,stsm-method}
\alias{set.sgfc}
\alias{set.sgfc,stsm-method}
\alias{set.xreg}
\alias{set.xreg,stsm-method}

\title{Setter Methods for Class \code{stsm}}

\description{
Setter or modifier methods for objects of class \code{\link{stsm}}.
}

\usage{
\S4method{set.cpar}{stsm}(x, value, check = TRUE, inplace = FALSE)
\S4method{set.nopars}{stsm}(x, v, check = TRUE, inplace = FALSE)
\S4method{set.pars}{stsm}(x, v, check = TRUE, inplace = FALSE)
\S4method{set.sgfc}{stsm}(x, inplace = FALSE)
\S4method{set.xreg}{stsm}(x, xreg, coefs = NULL)
}

\arguments{
\item{x}{an object of class \code{\link{stsm}}.}
\item{value}{a numeric value.}
\item{v}{a numeric vector.}
\item{check}{logical. If \code{TRUE}, the resulting model is checked for 
consistency with the definition of the \code{stsm} object.}
\item{inplace}{logical. If \code{TRUE}, the input object \code{x} is modified 
in place instead of returning the whole object.}
\item{xreg}{a matrix or numeric vector of external regressors. The number of 
rows or length of the vector must be equal to the length of \code{x@y}. 
If column names are specified they are used to name the parameters in 
the slot \code{pars}.}
\item{coefs}{an optional vector containing the value of the coefficients 
related to the regressors \code{xreg}. If the elements of the vector do not 
contain names they are assumed to be defined in the same order as the columns 
in the matrix \code{xreg}.}
}

\details{
\emph{Models parameterized with non-null \code{transPars}.}
If the model is parameterized according to a non-null value 
of the slot \code{transPars}, the argument \code{v} 
must contain the values of the auxiliary set of parameters \eqn{\theta}
rather than the actual parameters (variances and autoregressive coefficients).
For example, with \code{x@transPars = "square"} the variances are \eqn{\theta^2}.
Although this design may seem to disagree with the 
getter methods \code{\link{stsm-get-methods}}, 
the relevant input for the setter methods is actually the auxiliary values \eqn{\theta}. 
Be aware that if \code{transPars} is not null the parameters are 
transformed by \code{get.pars} according to the selected parameterization. 
Therefore, \code{v} must be referred to the non-transformed parameters.
%or the inverse of the transformed parameters

The previous comment does not apply to the argument \code{value} since \code{cpar} 
is not affected by \code{transPars}.
%It is assumed that suitable value is chosen, for example
%the value that maximizes the concentrated likelihood function.

\emph{Setter methods are safer.}
For those users that are not familiar with the design and internal structure of the 
class \code{\link{stsm}}, it is safer to use setter methods 
rather than modifying the contents of the slots through the \code{@<-} operator. 
See the examples below.
%\code{set.pars{x, v}} is a safer alternative to \code{x@pars <- v}.

\emph{Modifying the input object in-place.}
Instead of returning the whole object and create a new one or overwrite 
the original, it is possible to modify just the desired slot in the original
object that is passed as input. In the former case the \code{stsm} object returned by 
the method must be assigned to another object using the usual operator \code{<-}.
In the latter approach, the \code{stsm} object that is passed as argument is modified 
in-place. See the example below. The solution to modify an object in-place is taken from 
\href{http://tolstoy.newcastle.edu.au/R/help/04/02/0966.html}{this post}.
This option is not a customary solution in \R, however, it seems suitable in this context. 
The real benefit of this approach would depend on how \R deals with objects that 
are returned from functions. If assigning the output to a new object involves making 
copies of all the slots, then modifying the object in-place would most likely be more efficient 
since the desired slot is directly modified avoiding copying the whole object. 

After R version 3.1 this issue may become less critical. 
One of the new features reported in the release of R 3.1 states:
\emph{Avoid duplicating the right hand side values in complex assignments when possible. 
This reduces copying of replacement values in expressions such as Z$a <- a0}.
A related discussion for S4 classes can be found in 
\href{http://stackoverflow.com/questions/22448198/does-r-copy-unevaluated-slots-in-s4-classes-on-assignment/}{this post}.

\emph{Constant terms in the spectral generating function.}
In pure variance models, some elements of the spectral generating function (s.g.f.) 
do not depend on the parameters and can be stored as constants.
The method \code{set.sgfc} computes and stores those elements as a matrix
in the slot \code{sgfc}. This is useful for example when working with maximum likelihood methods 
in the frequency domain. In that context, the spectral generating function 
has to be updated several times for different parameter values. Having the information about 
the constant terms in the slot \code{sgfc} saves several computations whenever 
the s.g.f. is requested. For details about the s.g.f see \code{\link{stsm.sgf}}.

\emph{Further setter methods.}
Future versions may include additional setter methods, for example to 
change the slot \code{model} or to modify the time series \code{x@y}. 
The latter would also require updating the slots \code{diffy} and \code{ssd} 
if requested. Additional methods are not available in the current version because  
defining a new object by means of \code{\link{stsm.model}} will
often be better than modifying one of those slots that do not 
have a setter method.
}

\value{
If the slot is modified in place, \code{inplace=TRUE},
nothing is returned, the corresponding slot of the object \code{m} 
passed as argument is modified in place.

If \code{inplace=FALSE}, a new \code{stsm} object is returned.
It contains the same information as the input object \code{m} 
except for the slot that has been modified.
}

\seealso{
\code{\link{stsm}},
\code{\link{stsm.sgf}}.
}

\examples{
# sample models with arbitrary parameter values
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 15, "var3" = 30))
get.pars(m)

# correct modification
m1 <- set.pars(m, c(1, 2, 3))
get.pars(m1)
m1 <- set.pars(m, c(var1 = 11))
get.pars(m1)

# correct but error prone
m1@pars[] <- c(4, 22, 33)
get.pars(m1)
m1@pars <- c(var1 = 1, var2 = 2, var3 = 3)
get.pars(m1)

# inconsistent assignment (error returned)
# 'var4' is not a parameter of model 'llm+seas'
try(m1 <- set.pars(m, c(var4 = 4)))
# inconsistent assignment (no error returned)
# the error is not noticed at this point
# unless 'validObject' is called
m1 <- m
m1@pars["var4"] <- 4
get.pars(m1)
try(validObject(m1))

# modify only one element
m1 <- set.pars(m, v=c(var1=22))
get.pars(m1)
# wrong assignment, the whole vector in the slot is overwritten 
# no error returned at the time of doing the assignment
m1@pars <- c(var1 = 1)
get.pars(m1)
try(validObject(m1))

# consistent assignment but maybe not really intended
# all the elements are set equal to 12
m1 <- m
m1@pars[] <- 12
get.pars(m1)
# warning returned by 'set.pars'
m2 <- set.pars(m, 12)
get.pars(m2)

# wrong value unnoticed (negative variance)
m1 <- m
m1@pars[] <- c(-11, 22, 33)
get.pars(m1)
# negative sign detected by 'set.pars'
try(m1 <- set.pars(m, c(-11, 22, 33)))

# inplace = FALSE
# the whole object 'm' is assigned to a new object, 
# which will probably involve making a copy of all the slots
m <- set.pars(m, c(1,2,3), inplace = FALSE)
get.pars(m)

# inplace = TRUE
# the output is not assigned to a new object 
# the only operation is the modification of the slot 'pars'
# no apparent additional internal operations such as copying unmodified slots
get.pars(m)
set.pars(m, c(11,22,33), inplace = TRUE)
get.pars(m)

# set a matrix of regressors
xreg <- cbind(xreg1 = seq_len(84), xreg2 = c(rep(0, 40), rep(1, 44)))
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, xreg = xreg)
m
# set a new matrix of regressors to an existing
xreg3 <- seq(length(m@y))
m2 <- set.xreg(m, xreg3)
m2
# remove the external regressors
m3 <- set.xreg(m, NULL)
m3
m3@xreg
# initialize the coefficients to some values
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson,
  pars = c("xreg1" = 10), xreg = xreg)
m
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson,
  pars = c("xreg2" = 20, "xreg1" = 10), xreg = xreg)
m
}

\keyword{methods}
