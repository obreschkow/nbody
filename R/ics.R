ics.earth = function() {
  m = c(cste$Msun,cste$Mearth) # [m] masses of the sun and earth
  x = cbind(c(-cste$Mearth/sum(m)*cste$AU,0,0),
            c(cste$Msun/sum(m)*cste$AU,0,0)) # [m] position matrix
  v = cbind(c(0,2*pi*x[1,1]/cste$yr,0),
            c(0,2*pi*x[1,2]/cste$yr,0)) # [m/s] velocity matrix
  return(list(m=m, x=x, v=v))
}

ics.halley = function() {
  s = 17.834*cste$AU # [m] semi-major axis of Halley's Comet
  e = 0.96714 # [-] orbital eccentricity
  x.peri <<- (1-e)*s # [m] perihelion distance; as global variable for later use
  x.ap <<- (1+e)*s # [m] aphelion distance
  v.ap = sqrt((1-e)*cste$G*cste$Msun/(1+e)/s) # [m/s] aphelion velocity
  P <<- 2*pi*s*sqrt(s/cste$G/cste$Msun) # [s] orbital period
  m = c(cste$Msun,2.2e14) # [m] masses of the sun and the comet
  x = cbind(c(0,0,0),c(x.ap,0,0)) # [m] position matrix
  v = cbind(c(0,0,0),c(0,v.ap,0)) # [m/s] velocity matrix
  t = 0
  return(list(m=m, x=x, v=v, t=t))
}
