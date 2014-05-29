\name{sim-data}
\docType{data}
\alias{sim-data}
\alias{llmseas}

\title{Simulated Data}

\description{
Time series simulated from the local level plus seasonal model.
}

\usage{
llmseas
}

\format{
A time series of length \eqn{120}.
}

\source{
The scripts to generate this and other series used in the vignette 
are stored in the files \file{stsm/inst/datagen-llm.R} and 
\file{stsm/inst/datagen-llmseas.R}. 
The first series returned by the latter script is the series that is 
loaded by \code{data("llmseas")}.
}

\keyword{datasets}
