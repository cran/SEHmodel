% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecoToxic.R
\name{plotEcotoxic}
\alias{plotEcoToxic}
\alias{plotEcotoxic}
\title{Plot internal toxic concentration method}
\usage{
plotEcoToxic(objectL, objectI, objectT, numind = 8)
}
\arguments{
\item{objectL}{A Landscape object}

\item{objectI}{An Individuals object}

\item{objectT}{A ToxicIntensityRaster, a 3d array of Toxic intensity over the time [t,x,y], (first indice is time) see \code{\link{toxicIntensity}}}

\item{numind}{An individual ID}
}
\description{
Plot a time series of internal toxic concentration for a given individual.
}

