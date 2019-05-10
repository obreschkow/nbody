#' @title Physical constants
#'
#' @description Useful physical constants in SI units.
#'
#' @details Some of these constants are used internally by routines of the \code{nbody} package. Note that even if these constants are modified (e.g. via \code{cst$G=1}), the routines of the \code{nbody} package will still access the original values.
#'
#' @examples
#' summary(cst)
#'
#' @author Danail Obreschkow
#'
#' @export cst

cst = list(G = 6.67408e-11, # [m^3/kg/s^2] gravitational constant
           Msun = 1.98847e30, # [kg] solar mass
           Mearth = 5.972e24, # [kg] Earth mass
           AU = 149597870700, # [m] astronomical unit
           yr = 31557600, # [s] Julian astronomical year year
           year = 31557600, # [s] Julian astronomical year year
           month = 2629743.83333, # [s] month = year/12
           week = 7*86400, # [s] week
           day = 86400, # [s] day
           hour = 3600, # [s] hour
           min = 60, # [s] minute
           pc = 3.0857e16, # [m] parsec
           kpc = 3.0857e19, # [m] kiloparsec
           Mpc = 3.0857e22, # [m] Megaparsec
           Gpc = 3.0857e25, # [m] Gigaparsec
           lightyear = 9.4605284, # [m] light-year
           ly = 9.4605284) # [m] light-year
