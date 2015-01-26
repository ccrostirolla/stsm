\name{stsm-package}
\docType{package}
\alias{stsm-package}

\title{Structural Time Series Models}

\description{
This package provides algorithms to fit structural time series models 
by maximum likelihood.}

\details{
As witnessed in the special issue of the \emph{Journal of Statistical Software}
(Commandeur \emph{et al.}, 2011),
the prevalent procedure to fit a structural time series model is as follows:
1) choose arbitrary starting values for the parameters, 
2) evaluate the log-likelihood function by means of the Kalman filter, 
3) obtain a new set of parameter values that lead to a higher value
of the log-likelihood function by means of the L-BFGS-B algorithm,
4) iterate the searching procedure until a predetermined degree of convergence.
%
Considering that there are several packages in \R to run the Kalman filter 
(see for instance Tusell, 2011 and the documentation in package \pkg{KFKSDS})
and that the \code{\link[stats]{optim}} function in the \pkg{stats} package provides 
and interface to the L-BFGS-B and to other optimization algorithms,
fitting a structural time series model may seem a simple procedure that 
requires little more than translating the matrices of the state space form of 
the model into the syntax of the chosen interface.

In practice, the process is not always that straightforward.
As stated in the documentation of \code{\link[stats]{StructTS}},
\emph{optimization of structural models is a lot harder than many of the 
references admit}.
There are several details that should be taken into account 
when implementing the procedure described above, 
(López-de-Lacalle, 2013). 

There are not many packages in \R that provide alternative procedures to fit 
structural models. It is probably a consequence of the widespread believe that 
all that is needed to carry out and analysis with structural time series models 
is an implementation of the Kalman filter together with a general purpose
optimization algorithm.

%The package \pkg{stsm} provides algorithms to fit structural time series 
%models by maximum likelihood. An interface to the common procedure described above is 
%provided. It includes some choices that are often omitted or ignored. In addition, 
%other optimization procedures are provided: a scoring algorithm that maximizes 
%the frequency domain likelihood function and the Expectation-Maximization algorithm. 
%When possible, the algorithms are enhanced with analytical derivatives.
%Thus, the package \pkg{stsm} provides algorithms specialized to fit structural time 
%series models by maximum likelihood.
%Upon the utilities available in the packages \pkg{stsm.class} and \pkg{KFKSDS}, 
The package \pkg{stsm} implements specific algorithms to fit models in 
the framework of the basic structural time series model. 
The following enhancements to general purpose optimization algorithms are implemented:
scoring algorithm based on analytical derivatives, 
maximization of the time or frequency domain likelihood function,
automatic choice of the optimal step size, concentration of a parameter,
implementation of the original and a modified version of the 
expectation-maximization algorithm. 
}

\author{
Javier López-de-Lacalle \email{javlacalle@yahoo.es}
}

\section{References}{
Commandeur, J.J.F., Koopman, S.J. and Ooms, M. (2011).
\sQuote{Statistical Software for State Space Methods},
\emph{Journal of Statistical Software}, Vol. 41, No. 1,
\url{http://www.jstatsoft.org/v41/i01/}.

Durbin, J. and Koopman, S.J. (2001).
\emph{Time Series Analysis by State Space Methods}. 
Oxford University Press.

Harvey, A.C. (1989).
\emph{Forecasting, Structural Time Series Models and the Kalman Filter}. 
Cambridge University Press.

López-de-Lacalle, J. (2013).
\sQuote{101 Variations on a Maximum Likelihood Procedure for a Structural Time Series Model.}
Unpublished manuscript.

Tusell. F. (2011).
\sQuote{Kalman Filtering in \R.} 
\emph{Journal of Statistical Software}, Vol. 39, No. 2. 
\url{http://www.jstatsoft.org/v39/i02/}.
}

\keyword{package, ts}
