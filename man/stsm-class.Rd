\name{stsm-class}
\docType{class}
\alias{stsm-class}
\alias{stsm}

\title{Class \code{stsm} for Structural Time Series Models}

\description{
This class defines a structural time series model.
}

\section{Slots}{
\describe{
%
\item{\code{call}}{Object of class \code{call}. Call to \code{\link{stsm.model}}.}
%
\item{\code{model}}{Object of class \code{character}. Name or label for the selected model (see
\code{\link{stsm.model}} for available models).}
%
\item{\code{y}}{Object of class \code{ts}. Original time series.}
%
\item{\code{diffy}}{Object of class \code{ts}. Differenced series \code{y}. The differencing
operator that renders stationarity in the \code{model} is applied to the series \code{y}.}
%
\item{\code{xreg}}{An optinal matrix or numeric vector of external regressors.}
%
\item{\code{fdiff}}{Object of class \code{function}. Function with arguments 
\code{x}: a \code{ts} object, \code{s}: periodicity of the data. This function
applies the differencing operator that renders stationarity in the \code{model} 
to a \code{ts} object passed to it.}
%
\item{\code{ss}}{Object of class \code{list}. Matrices of the state space form 
of the structural model.}
%
\item{\code{pars}}{Object of class \code{numeric}. Named vector with the parameters of the model.}
%
\item{\code{nopars}}{Optional object of class \code{numeric}. An optional named vector with the remaining 
parameters of the model not included in \code{pars}.
This slot is not affected by the transformation of parameters \code{\link{transPars}}.
These parameters are considered fixed in the optimization procedures implemented in package \pkg{stsm}.}
%
\item{\code{cpar}}{Optional object of class \code{numeric}. Named vector of length one containing the
parameter that is concentrated out of the likelihood function (if any).}
%
\item{\code{lower}}{Object of class \code{numeric}. Named vector with the lower bounds for 
\code{pars}.}
%
\item{\code{upper}}{Object of class \code{numeric}. Named vector with the upper bounds 
for \code{pars}.}
%
\item{\code{transPars}}{Character string referring to the parameterization of the model, see
\code{\link{transPars}}.}
%
\item{\code{ssd}}{Optional object of class \code{numeric}. Sample spectral density 
(periodogram) of the differenced series \code{diffy}.}
%
\item{\code{sgfc}}{Optional object of class \code{matrix}. Constant elements in the spectral 
generating function of the \code{model} (for pure variance models).}
}}

\section{Methods}{
\describe{
%
\item{\code{char2numeric}}{Return a list containing the matrices of the state
space representation of the model. The matrices are the same as those in the slot \code{ss}
but the characters are replaced by the corresponding numeric values defined 
in \code{pars}, \code{nopars} and \code{cpar}.}
%
\item{\code{checkbounds}}{Check whether the values of \code{pars} lie within the 
\code{lower} and \code{upper} bounds.}
%
\item{\code{get.pars}}{Return the slot \code{pars}, the parameters of the model.
If the model is parameterized in terms of a set of auxiliary parameters such as those considered 
in \code{\link{transPars}}, then the transformed parameters are returned.
Thus, when the slot \code{transPars} is not \code{NULL} 
\code{x@pars} will not be equal to \code{get.pars(x)}.}
%
\item{\code{get.cpar}}{Return the slot \code{cpar}.}
%
\item{\code{get.nopars}}{Return the slot \code{nopars}.}
%
\item{\code{set.cpar}}{Set or modify the value of the slot \code{cpar}}
%
\item{\code{set.nopars}}{Set or modify the value of the slot \code{nopars}.}
%
\item{\code{set.pars}}{Set or modify the value of the slot \code{pars}.}
%
\item{\code{set.sgfc}}{Compute and set the value of the slot \code{sgfc}.}
%
\item{\code{set.xreg}}{Set or modify the value of the slot \code{xreg}.}
%
\item{\code{setValidity}}{Check the validity of the arguments passed to the function.}
%
\item{\code{show}}{Show a brief summary of the object.}
%
\item{\code{transPars}}{Transform the parameters of the model according to the parameterization 
defined in the slot \code{transPars}.}
}}

\seealso{
\code{\link{stsm.model}}.
}

\section{References}{
Christophe Genolini. \emph{A (Not So) Short Introduction to S4.
Object Oriented Programming in \R.} V0.5.1. August 20, 2008.
}

\keyword{classes}
