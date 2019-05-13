% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reset.cm.R
\name{reset.cm}
\alias{reset.cm}
\title{Move center of mass to the origin}
\usage{
reset.cm(sim)
}
\arguments{
\item{sim}{list of m, x, v or list with a sublist "ics", made of m, x, v, where
\code{m} = N-vector with the masses of the N particles\cr
\code{x} = N-by-3 matrix specifying the initial position in cartesian coordinates\cr
\code{v} = N-by-3 matrix specifying the initial velocities\cr\cr}
}
\value{
Returns a structure of the same format as the input argument, but with recentered positions and velocities.
}
\description{
Move center of mass to the origin and cancel center of mass velocity
}
\author{
Danail Obreschkow
}