#' @title Gravitational N-body simulation
#'
#' @description This package contains tools to run simple direct gravitational N-body simulations. It uses a variable block timestep and lets the user choose between a range of integrators, including 4th and 6th order integrators for high-accuracy simulations. Basic top-hat smoothing is available as an option. The code also allows the definition of background particles that are fixed or in uniform motion, not subject to acceleration by other particles.
#'
#' @author Danail Obreschkow <danail.obreschkow@icrar.org>
#'
#' @docType package
#' @name nbody
#' @useDynLib nbody
NULL
