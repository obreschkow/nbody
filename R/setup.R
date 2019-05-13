#' @title Initialize N-body simulation
#'
#' @description Routines to generate the structured lists of initial conditions and simulation parameters required to run an N-body simulation with \code{\link{run.simulation}}.
#'
#' @name setup
#' @aliases setup.earth
#' @aliases setup.halley
#'
#' @param t.max final simulation time
#' @param dt.out time step for simulation output
#' @param nperiods number of orbital periods to be computed (only for \code{setup.halley}); ignored if \code{t.max} is specified.
#' @param e eccentricity (only for \code{setup.halley})
#' @param s semi-major axis (only for \code{setup.halley})
#' @param ... other simulation parameters used by \code{\link{run.simulation}}
#'
#' @examples
#' sim = setup.halley()
#' sim = run.simulation(sim)
#' plot(sim)
#'
#' @author Danail Obreschkow
#'
#' @rdname setup
#' @return default routine, which is identical to setup.halley()
#' @export
setup <- function() setup.halley()

#' @rdname setup
#' @return \code{setup.sunearth()} sets up a simple 2-body simulation of the Earth around the Sun, using only approximate orbital specifications.
#' @export
setup.sunearth = function(t.max=cst$yr, dt.out=cst$yr, ...) {

  m = c(cst$Msun,cst$Mearth) # [m] masses of the sun and earth
  x = rbind(c(-cst$Mearth/sum(m)*cst$AU,0,0),
            c(cst$Msun/sum(m)*cst$AU,0,0)) # [m] position matrix
  v = rbind(c(0,2*pi*x[1,1]/cst$yr,0),
            c(0,2*pi*x[2,1]/cst$yr,0)) # [m/s] velocity matrix

  sim = list(ics = list(m=m, x=x, v=v),
             para = list(t.max = cst$yr, dt.out = cst$week, ...))

  class(sim) = 'simulation'
  return(sim)

}

#' @rdname setup
#' @return \code{setup.halley()} sets up a 2-body simulation of Halley's Commet around the Sun.
#' @export
setup.halley = function(t.max=NULL, nperiods=1, dt.out=3*cst$month, e=0.96714, s=17.834*cst$AU, ...) {

  s = 17.834*cst$AU # [m] semi-major axis of Halley's Comet
  e = 0.96714 # [-] orbital eccentricity
  x.peri = (1-e)*s # [m] perihelion distance; as global variable for later use
  x.ap = (1+e)*s # [m] aphelion distance
  v.ap = sqrt((1-e)*cst$G*cst$Msun/(1+e)/s) # [m/s] aphelion velocity
  period = 2*pi*s*sqrt(s/cst$G/cst$Msun) # [s] orbital period
  m = c(cst$Msun,2.2e14) # [m] masses of the sun and the comet
  x = rbind(c(0,0,0),c(x.ap,0,0)) # [m] position matrix
  v = rbind(c(0,0,0),c(0,v.ap,0)) # [m/s] velocity matrix

  sim = list(ics = list(m=m, x=x, v=v),
             para = list(t.max = nperiods*period, dt.out = dt.out, ...),
             user = list(period = period, x.peri = x.peri, x.ap = x.ap))

  class(sim) = 'simulation'
  return(sim)
}
