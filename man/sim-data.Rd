\name{sim-data}
\docType{data}
\alias{sim-data}
\alias{llmseas}
\alias{llm}

\title{Simulated Data Sets}

\description{
Time series simulated from the local level model and the local level plus seasonal model.
Each data set consists of \eqn{1000} series of \eqn{120} observations.
These data are used in the package's vignette to illustrate the performance of 
the algorithms.
}

\usage{
llm
llmseas
}

\format{
A \eqn{1000\times 120} matrix containing by column each simulated series.
}

\source{
The scripts used to generate the data are stored in the files
\file{stsm/inst/datagen-llm.R} and \file{stsm/inst/datagen-llmseas.R}.
The series in the data set \code{llmseas} are simulated by means of 
\code{\link{datagen.stsm}}. See the examples in the documentation of 
that function for details.
}

\keyword{datasets}
