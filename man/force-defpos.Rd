\name{force.defpos}
\alias{force.defpos}

\title{Force Positive Definiteness of a Matrix}

\description{
This function transforms a matrix to be positive definite.
}

\usage{
force.defpos(m, tol = 0.001, debug = FALSE)
}

\arguments{
\item{m}{a matrix.}
\item{tol}{tolerance.}
\item{debug}{logical. If \code{TRUE}, it is checked that the output matrix 
is positive definite.}
}

\details{
%A negative definite matrix is required in the Newton-Raphson or the scoring algorithms. 
%It ensures that the gradient is projected direction towards the location of a local maximum.
%Despite in theory the Hessian is definite positive, in practice it may turn to be 
%the opposite.
%
The scoring algorithm \code{\link{maxlik.fd.scoring}} requires a positive definite matrix
to project the gradient into the optimal direction. 
If that matrix happens to be non-positive definite then the matrix 
\eqn{M} is transformed as described in Pollock (1999) pp. 341-342:
%
\deqn{M = M + (\mu - \kappa) I}
%
where \eqn{I} is the identity matrix, \eqn{\mu} is a tolerance value and 
\eqn{\kappa} is the smallest eigenvalue of \eqn{M}.
%
Other alternatives are discussed in Nocedal and Wright (2006) chapter 3.
}

\value{
A positive definite matrix.
}

\references{
Nocedal, J. and Wright, J. W. (2006).
\emph{Numerical Optimization}. Springer-Verlag.

Pollock, D.S.G. (1999).
\emph{A Handbook of Time-Series Analysis Signal Processing and Dynamics}.
Academic Press.
}

\keyword{algebra}
\keyword{array}
