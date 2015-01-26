\name{stsm.model}
\alias{stsm.model}

\title{Wrapper for Constructor of Objects of Class \code{stsm}}

\description{
Interface to define an object of class \code{\link{stsm}}.
This is a wrapper function to constructor \code{\link[methods]{new}}.
}

\usage{
stsm.model(model = c("local-level", "local-trend", "BSM", 
  "llm+seas", "trend+ar2"),
  y, pars = NULL, nopars = NULL, cpar = NULL, xreg = NULL, 
  lower = NULL, upper = NULL, transPars = NULL,
  ssd = FALSE, sgfc = FALSE)
}

\arguments{
\item{model}{a character selecting the structural time series model.}
%
\item{y}{a univariate time series, \code{\link[stats]{ts}}.}
%
\item{pars}{initial values for the parameters of the model. It must be a named vector.}
%
\item{nopars}{optional named numeric containing the remaining parameters of the model 
not included in \code{pars} and \code{cpar}.}
%
\item{cpar}{optional named numeric of length one. See details.}
%
\item{xreg}{optional matrix or numeric vector of external regressors.}
%
\item{lower}{optional named vector setting lower bounds to some parameters of the model. 
The names must follow the same same labelling as \code{pars}.}
%
\item{upper}{optional named vector setting upper bounds to some parameters of the model. 
The names must follow the same same labelling as \code{pars}.}
%
\item{transPars}{optional character choosing one of the parameterizations
defined in \code{\link{transPars}} or a function defining an alternative
parameterization.}
%
\item{ssd}{logical. If \code{TRUE}, the sample spectral density (periodogram) 
of the stationary transformation of the data is computed and stored in the slot \code{ssd}.
Otherwise, it is ignored.}
%
\item{sgfc}{logical. If \code{TRUE}, constants terms of the spectral generating
function related to the chosen model are computed and stored in the slot \code{sgfc}.
Otherwise, it is ignored.}
}

\section{Available models}{
The \bold{local level model} consists of a random walk plus a Gaussian disturbance term.

The measurement equation is:
\deqn{y[t] = m[t] + e[t], e[t] \sim  N(0, \sigma^2_1)}

The state equation is:
\deqn{m[t+1] = m[t] + v[t], v[t] \sim N(0, \sigma^2_2)}

The \bold{local trend model} consists of a trend where the slope evolves as
a random walk.

The measurement equation is:
\deqn{y[t] = m[t] + e[t], e[t] \sim  N(0, \sigma^2_1)}

The state equations are:
\deqn{m[t+1] = m[t] + n[t] + v[t], v[t] \sim N(0, \sigma^2_2)}
\deqn{n[t+1] = n[t] + w[t], w[t] \sim N(0, \sigma^2_3)}

Setting \eqn{var3 = 0} yields the local level model. 
The constraint \eqn{var2 = 0} involves a smooth trend.

The \bold{basic structural model} consists of a local trend model 
plus a seasonal component. 

The measurement equation is:
\deqn{y[t] = m[t] + s[t] + e[t], e[t] \sim N(0, \sigma^2_1)}

The state equations are the same as the local trend model plus
a seasonal component:
\deqn{s[t+1] = -s[t] - \dots - s[t-freq+2] + w[t],
w[t] \sim N(0, \sigma^2_4)}

The restriction \eqn{\sigma^2_4 = 0} yields a deterministic seasonal pattern.

According to the labelling convention used in the package,
the variance parameters \eqn{\sigma^2_1}, \eqn{\sigma^2_2}, \eqn{\sigma^2_3}
and \eqn{\sigma^2_4} are respectively denoted
\code{"var1"}, \code{"var2"}, \code{"var3"} and \code{"var4"}.
}

\details{
\emph{Slot \code{pars} and \code{nopars}.}
In some situations it is convenient to split the vector of parameters in 
two vectors, the slot \code{pars} and the slot \code{nopars}. For example, 
when the parameters are to be estimated by an optimization algorithm, 
only the parameters in \code{pars} are allowed to change while the parameters in
\code{nopars} are considered fixed.

\emph{Scaling parameter \code{cpar}.}
The model can be defined in terms of relative variances.
In this case, the variance that acts as a scaling parameter is stored in the slot \code{cpar}. 
Otherwise, \code{cpar} is null and ignored.
Typically, the scaling parameter will be chosen to be the variance parameter that is 
concentrated out of the likelihood function.

\emph{Naming convention of parameters.}
The parameters defined in the slots \code{pars}, \code{nopars} and \code{cpar} 
must be labelled according to the following naming convention. 
The variance parameters abide by the regular expression
\dQuote{^var\\d\{1,2\}$}, e.g. \code{var1}, \code{var2},... 
The variances of the initial state vector, \eqn{P0}, follow a similar 
naming convention, \code{P01}, \code{P02},...
The elements of the initial state vector, \eqn{a0}, are similarly denoted
as \code{a01}, \code{a02},...

\emph{Default values.}
Default values are assigned to the slots \code{pars}, \code{nopars} and \code{cpar} if they 
are no defined in their corresponding arguments passed to \code{stsm.model}.
By default, the variance parameters are defined in the slot \code{pars} with value \eqn{1}.
The initial state vector is assigned by default to \code{nopars}, it takes on zero values 
except for the first element that takes the value of the 
first observation in the data. The variance of the initial state vector is assigned 
by default to \code{nopars} as well. By default it takes on the value \eqn{10000} 
times the variance of the data.

If the argument \code{pars} is not \code{NULL}, no other parameters are stored in 
the slot \code{pars}. If the argument \code{nopars} is not \code{NULL}, the parameters
in that argument are added to the other default parameters. 
This is more convenient in practice. See the examples below.

\emph{Alternative parameterizations.}
See \code{\link{transPars}} for available parameterizations of 
the model. The definition of a function to be defined in the slot \code{transPars} is 
also explained there.

\emph{Stationary transformation of the data.}
The sample spectral density is computed for the differenced time series \code{y}.
The differencing filter is chosen so that the data are rendered stationary according 
to the selected model. The stationary form of each model is given in 
\code{\link{stsm.sgf}}.
}

\value{
An object of class \code{stsm}.
}

\seealso{
\code{\link{stsm}}.
}

\examples{
# sample model with arbitrary parameter values
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c("var1" = 2, "var2" = 6), nopars = c("var3" = 12))
m

# parameter values
v <- c("var1" = 2, "var2" = 6, "var3" = 3, "var4" = 12)

# define the parameter 'cpar' and let the the remaining parameters 
# be defined by default in the slots 'pars' and 'nopars'
m <- stsm.model(model = "BSM", y = JohnsonJohnson, 
  pars = NULL, nopars = NULL, cpar = v[1])
m@pars
m@nopars
m@cpar

# define the slot 'pars', only 'v[1]' is stored in 'pars'
# the remaining variances are moved to 'nopars' along
# with the initial state vector and its variances
m <- stsm.model(model = "BSM", y = JohnsonJohnson, 
  pars = v[1])
m@pars
m@nopars
m@cpar

# define some of the parameters to be stored in the slot 'nopars' 
# 'only 'v[1:2]' is added to the remaining elements in 'nopars'; 
# by default the variances not defined in 'nopars' are assigned to 
# 'pars' with value 1
m <- stsm.model(model = "BSM", y = JohnsonJohnson, 
  nopars = v[1:2])
m@pars
m@nopars
m@cpar

# define the slot 'pars' and set a particular value to 
# some variances stored in 'nopars', 'v[2:3]'
# 'var4' takes the default value 1 and is stored in 'nopars' 
# since the definition 'pars = v[1]' excludes it form 'pars'
m <- stsm.model(model = "BSM", y = JohnsonJohnson, 
  pars = v[1], nopars = v[2:3])
m@pars
m@nopars
m@cpar

# define the slots 'pars' and 'cpar'
# the remaining parameters are stored in 'nopars' with the 
# values by default
m <- stsm.model(model = "BSM", y = JohnsonJohnson, 
  pars = v[2:4], nopars = NULL, cpar = v[1])
m@pars
m@nopars
m@cpar
}

\keyword{ts, model}
