% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nbody-package.R
\docType{package}
\name{nbody-package}
\alias{nbody}
\alias{nbody-package}
\title{Gravitational N-Body Simulation}
\description{
Run simple direct gravitational N-body simulations. The package can access different external N-body simulators (e.g. GADGET-4 by Springel et al., 2021), but also has a simple built-in simulator. This default simulator uses a variable block time step and lets the user choose between a range of integrators, including 4th and 6th order integrators for high-accuracy simulations. Basic top-hat smoothing is available as an option. The code also allows the definition of background particles that are fixed or in uniform motion, not subject to acceleration by other particles.
}
\author{
Danail Obreschkow <danail.obreschkow@gmail.com>
}
