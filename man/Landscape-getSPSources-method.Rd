% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Landscape-Methods.R
\docType{methods}
\name{getSPSources}
\alias{Landscape,getSPSources-method}
\alias{getSPSources}
\alias{getSPSources,Landscape-method}
\title{Method to get SpatialPolygons Sources from a Landscape}
\usage{
getSPSources(object, ...)

\S4method{getSPSources}{Landscape}(object)
}
\arguments{
\item{object}{A Landscape object}

\item{...}{options}
}
\value{
a SpatialPolygonDataFrame object
}
\description{
Method to get SpatialPolygons Sources from a Landscape

Get all polygons of a Landscape object identified as sources
}
\examples{
\dontrun{
## To get the first spatial polygon source information from a landscape object (named land)
getSPSources(land)@polygons[[1]]
}
}

