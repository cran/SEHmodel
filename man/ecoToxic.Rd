% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecoToxic.R
\docType{methods}
\name{ecoToxic}
\alias{ecoToxic}
\alias{ecoToxic,Landscape,Individuals-method}
\alias{ecotoxicological-method}
\title{EcoToxicological model method}
\usage{
ecoToxic(objectL, objectI, ...)

\S4method{ecoToxic}{Landscape,Individuals}(objectL, objectI, objectT,
  mintime = 1, maxtime, kin = 0.25, kout = 0.5, deltat = 0.1)
}
\arguments{
\item{objectL}{A Landscape object}

\item{objectI}{An Individuals object}

\item{...}{other parameters}

\item{objectT}{A ToxicIntensityRaster object, a 3D array of Toxic Dispersion over time [t,x,y], first indice is time}

\item{mintime}{Time to start simulation (default = 1)}

\item{maxtime}{Time to end simulation}

\item{kin}{ingestion rate (\% of contaminants staying in the body)}

\item{kout}{elimination rate (\% of contaminants eliminated from the body)}

\item{deltat}{\% of a time unit for the ordinary differential equation (ODE)}
}
\value{
An Individuals object with updated internal toxic concentrations
}
\description{
Generic method on Landscape and Individuals objects applying ecotoxicological equation.

This method gives internal concentration of contaminants within individuals, from toxic quantities in the environment and individual parameters.
}

