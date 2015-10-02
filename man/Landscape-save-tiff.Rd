% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Landscape-Methods.R
\name{saveIntoFIle}
\alias{saveIntoFIle}
\alias{saveIntoFile}
\title{Save Particles Dispersion 3d Array to tiff file}
\usage{
saveIntoFile(objectL, objectT, filename = "ParticlesDispersion.tif",
  format = "GTiff", overwrite = T)
}
\arguments{
\item{objectL}{a Landscape object}

\item{objectT}{a 3d array particles dispersion indexed by time}

\item{filename}{output file name}

\item{format}{output format (default=GTiff)}

\item{overwrite}{if True overwrite filename}
}
\value{
RasterStack
}
\description{
Save into tiff file particles dispersion 3d array from toxicIntensity. Output a RasterStack with a layer by time step with projection set to CRS="+proj=longlat +datum=WGS84"
}
