\name{stsm-char2numeric-methods}
\docType{methods}
\alias{char2numeric}
\alias{char2numeric,stsm-method}

\title{State Space Representation of Objects of Class \code{stsm}}

\description{
This method returns the state space representation of time series models 
defined in the class \code{\link{stsm}}.
}

\usage{
\S4method{char2numeric}{stsm}(x, P0cov = FALSE, rescale = FALSE)
}

\arguments{
\item{x}{an object of class \code{\link{stsm}}.}
\item{P0cov}{logical. If \code{TRUE} the values of the elements outside the diagonal 
in the initial covariance matrix of the state vector are set equal to the values in the
diagonal. Otherwise values outside the diagonal are set equal to zero.}
\item{rescale}{logical. If \code{TRUE}, relative variance parameters are rescaled into 
absolute variances. Otherwise, relative variances are used. 
Ignored if \code{x@cpar} is null.}
}

\section{State space representation}{
The general univariate linear Gaussian state space model is defined 
as follows:
\deqn{y[t] = Za[t] + e[t], e[t] \sim~ N(0, H)}
\deqn{a[t+1] = Ta[t] + Rw[t], w[t] \sim~ N(0, V)}

for \eqn{t=1,\dots,n} and \eqn{a[1] \sim~ N(a0, P0)}.
\eqn{Z} is a matrix of dimension \eqn{1\times m}{1xm}; 
\eqn{H} is \eqn{1\times 1}{1x1};
\eqn{T} is \eqn{m\times m}{mxm};
\eqn{R} is \eqn{m\times r}{mxr};
\eqn{V} is \eqn{r\times r}{rxr}; 
\eqn{a0} is \eqn{m\times 1}{mx1} and
\eqn{P0} is \eqn{m\times m}{mxm},
where \eqn{m} is the dimension of the state vector \eqn{a} and 
\eqn{r} is the number of variance parameters in the state vector.
}

\details{
%See \code{\link{stsm.model}} for the notation used 
%in the state space representation.
%
This method uses the information from the slots \code{pars}, \code{nopars} 
and \code{cpar} in order to build the numeric representation of the matrices.

For details about the argument \code{rescale} see the details section 
in \code{\link{stsm-get-methods}} and the examples below.

A previous version of this method employed the information in the slot 
\code{ss}. This slot contains the matrices of the state space form of the model 
but instead of inserting the parameter values, character strings indicating 
the location of the parameters are placed in the corresponding cells. 
This method performed the mapping from the character to the numeric matrices 
by means of a internal function called \code{ss.fill}.
Currently the slot \code{ss} and the matrices are directly built depending on 
the model that was selected among those available in \code{\link{stsm.model}}.
The current approach is straightforward and faster. 
The previous approach may still be interesting to allow the user to define 
additional models just by translating the notation of the model into character
matrices. The usefulness of enhancing this approach will be assessed in 
future versions of the package.
}

\value{
A list of class \code{stsmSS} containing the following numeric matrices and vectors:
\item{Z}{observation matrix.}
\item{T}{transition matrix.}
\item{H}{observation variance.}
\item{R}{selection matrix.}
\item{V}{state vector variance-covariance matrix.}
\item{Q}{RVR'.}
\item{a0}{initial state vector.}
\item{P0}{initial state vector uncertainty matrix.}

The list contains also two vectors, \code{Vid} and \code{Qid}, 
with the indices of those cells where the variance parameters 
are located respectively in the matrices \eqn{V} and \eqn{Q}.
The first element in a matrix is indexed as \eqn{0}.
}

\seealso{
\code{\link{stsm}},
\code{\link{stsm.model}}.
}

\examples{
# sample model with arbitrary parameter values
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 6), nopars = c("var3" = 12))
ss1 <- char2numeric(m)
c(get.pars(m), get.nopars(m), get.cpar(m))
# character notation of the covariance matrix of the state vector
m@ss$Q
# information from the slots 'pars', 'nopars' and 'cpar'
# is used to retrieve the numeric representation of 'm@ss$Q'
ss1$Q

# same as above but with P0cov=TRUE
# the only change is in the initial covariance matrix of
# the state vector 'P0'
ss2 <- char2numeric(m, P0cov = TRUE)
ss1$P0
ss2$P0

# if a non-standard parameterization is used, 
# the values in the slot 'pars' are transformed accordingly 
# and the actual variance parameters are returned;
# notice that the transformation of parameters applies only 
# to the parameters defined in the slot 'pars'
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 6), nopars = c("var3" = 12),
  transPars = "square")
c(get.pars(m), get.nopars(m), get.cpar(m))[1:3]
ss <- char2numeric(m)
ss$H
ss$Q

# model defined in terms of relative variances,
# the variances in 'pars' are relative to the scaling parameter 'cpar',
# in this example 'cpar' is chosen to be the variance 'var1'
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var2" = 3, "var3" = 6), cpar = c("var1" = 2),
  transPars = NULL)
# the state space representation can be done with 
# relative variances (no rescaling)
ss <- char2numeric(m, rescale = FALSE)
ss$H
ss$Q
# or with absolute variances (rescaling)
ss <- char2numeric(m, rescale = TRUE)
ss$H
ss$Q

# in a model where the parameters are the relative variances 
# and with non-null 'transPars', the transformation is applied to 
# the relative variances, not to the absolute variances, i.e., 
# the relative variances are first transformed and afterwards they are 
# rescaled back to absolute variances if requested
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var2" = 3, "var3" = 6), cpar = c("var1" = 2),
  transPars = "square")
# the state space representation can be done with 
# relative variances (no rescaling)
ss <- char2numeric(m, rescale = FALSE)
ss$H
ss$Q
# or with absolute variances (rescaling)
ss <- char2numeric(m, rescale = TRUE)
ss$H
ss$Q
}

\keyword{methods}
