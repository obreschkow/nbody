#' Make ICs for N-body simulation
#'
#' Description
#'
#' @importFrom magicaxis magplot
#'
#' @param sim is an output of \code{run.simulation}
#'
#' @author Danail Obreschkow
#'
#' @export ics.earth ics.halley

ics.earth = function() {
  m = c(cst$Msun,cst$Mearth) # [m] masses of the sun and earth
  x = cbind(c(-cst$Mearth/sum(m)*cst$AU,0,0),
            c(cst$Msun/sum(m)*cst$AU,0,0)) # [m] position matrix
  v = cbind(c(0,2*pi*x[1,1]/cst$yr,0),
            c(0,2*pi*x[1,2]/cst$yr,0)) # [m/s] velocity matrix
  return(list(m=m, x=x, v=v))
}

ics.halley = function() {
  s = 17.834*cst$AU # [m] semi-major axis of Halley's Comet
  e = 0.96714 # [-] orbital eccentricity
  x.peri = (1-e)*s # [m] perihelion distance; as global variable for later use
  x.ap = (1+e)*s # [m] aphelion distance
  v.ap = sqrt((1-e)*cst$G*cst$Msun/(1+e)/s) # [m/s] aphelion velocity
  P = 2*pi*s*sqrt(s/cst$G/cst$Msun) # [s] orbital period
  m = c(cst$Msun,2.2e14) # [m] masses of the sun and the comet
  x = cbind(c(0,0,0),c(x.ap,0,0)) # [m] position matrix
  v = cbind(c(0,0,0),c(0,v.ap,0)) # [m/s] velocity matrix
  t = 0
  return(list(m=m, x=x, v=v, t=t,
              attributes = list(P = P, x.peri = x.peri, x.ap = x.ap)))
}
