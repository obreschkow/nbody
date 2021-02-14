#' Interpolate particle coordinates between snapshots
#' 
#' @importFrom graphics curve points
#' @importFrom plotrix draw.arc
#'
#' @description Interpolate/extrapolate particle positions and velocities between discrete time steps. An independent polynomial interpolation along each Cartesian coordinate is used by default. If an acceleration field is specified via the optional argument \code{afield}, a leapfrog integrator with the custom time step \code{dt} is applied. This is particularly useful for interpolation between time steps separated by one or more characteristic orbital times.
#' 
#' @param x0 n-by-d array of n particle positions in d dimensions at time t0; or an n-vector for one dimension
#' @param x1 n-by-d array of n particle positions in d dimensions at time t1; or an n-vector for one dimension
#' @param t0 single number giving the time of the first time step, corresponding to the particle positions x0
#' @param t1 single number giving the time of the second time step, corresponding to the particle positions x1
#' @param ti single number giving the time of the interpolated/extrapolated output positions and velocities
#' @param v0 optional n-by-d array of n particle velocities in d dimensions at time t0; or an n-vector for one dimension
#' @param v1 optional n-by-d array of n particle velocities in d dimensions at time t1; or an n-vector for one dimension
#' @param afield optional acceleration field for leapfrog integration. This must be a function of a n-by-d array representing n position vectors, which returns an n-by-d array with the n acceleration vectors in the same length and time units as x and v.
#' @param dt optional time step used for leapfrog integration; only used if \code{afield} given. The default is \code{dt=abs(t1-t0)/20}.
#'
#' @author Danail Obreschkow
#' 
#' @examples
#' 
#' ## Example 1: 1D-harmonic oscillator
#' 
#' # analytical model of continuous orbit in a 1D harmonic oscillator
#' x = function(t) sin(t) # position
#' v = function(t) cos(t) # velocity
#' curve(x,0,pi,200,xlab='Time t',ylab='Position x(t)')
#' 
#' # evaluation at two discrete times
#' t0 = 0.8; t1 = 2.8
#' points(c(t0,t1),x(c(t0,t1)),pch=16)
#' 
#' # interpolation/extrapolation
#' for (ti in seq(0,5,by=0.1)) {
#'   
#'   # using only positions (blue)
#'   p = particleinterp(x(t0),x(t1),t0,t1,ti)
#'   points(ti,p$x,col='blue')
#'   
#'   # using positions and velocities (red)
#'   p = particleinterp(x(t0),x(t1),t0,t1,ti,v(t0),v(t1))
#'   points(ti,p$x,col='red')
#'   
#' }
#' 
#' ## Example 2: particles moving in a spherical gravitational potential
#' 
#' # make acceleration field of an isothermal potential of unit circular velocity
#' afield = function(x) -x/rowSums(x^2)
#' 
#' # analytic function of circular orbits
#' x = function(t) cbind(x=r*cos(t/r+phi),y=r*sin(t/r+phi)) # positions
#' v = function(t) cbind(-sin(t/r+phi),cos(t/r+phi)) # velocities 
#' 
#' # specification of two circular orbits
#' r = c(1,2) # radii of circular orbits
#' phi = c(0,2) # azimuths at t=0 in radians
#' 
#' # times of the two time steps
#' t0 = 0 # initial time
#' t1 = 6 # final time
#' 
#' # plot initial and final positions as black points
#' plot(rbind(x(t0),x(t1)),pch=16,asp=1,xlim=c(-3,3),ylim=c(-3,3))
#' 
#' # draw lines for analytical circular orbits
#' plotrix::draw.arc(0,0,r,t0/r+phi,t1/r+phi)
#' 
#' # plot leapfrog-interpolated points as red crosses
#' for (ti in seq(t0,t1,length=20)) {
#'   out = particleinterp(x(t0),x(t1),t0,t1,ti,v(t0),v(t1),afield)
#'   points(out$x,col='red',pch=4)
#' }
#' 
#' # plot standard polynomial interpolation as blue circles
#' for (ti in seq(t0,t1,length=20)) {
#'   out = particleinterp(x(t0),x(t1),t0,t1,ti,v(t0),v(t1))
#'   points(out$x,col='blue')
#' }
#'
#' @export

particleinterp = function(x0,x1,t0,t1,ti=0.5,v0=NULL,v1=NULL,afield=NULL,dt=NULL) {
  
  st = ti-t0
  ht = t1-t0
  
  if (is.null(afield)) {
  
    if (is.null(v0)) v0 = (x1-x0)/ht
    if (is.null(v1)) v1 = (x1-x0)/ht
    
    if (ti<=t0) {
      
      # linear extrapolation before t0
      x = x0+v0*st
      v = v0
      
    } else if (ti>=t1) {
      
      # linear extrapolation after t1
      x = x1+v1*(ti-t1)
      v = v1
      
    } else {
      
      # polynomial interpolation between t0 and t1
      q = 2*(x0-x1)+(v0+v1)*ht
      p = 3*(x1-x0)-(2*v0+v1)*ht
      x = x0+v0*st+p*(st/ht)^2+q*(st/ht)^3
      v = v0+2*p*st/ht^2+3*q*st^2/ht^3
      
    }
    
  } else {
    
    if (is.null(v0)) stop('v0 must be given if afield is given')
    if (is.null(v1)) stop('v1 must be given if afield is given')
    
    if (is.null(dt)) dt = abs(t1-t0)/20
    
    leapfrog = function(x0,v0,tmax) {
      x = x0
      v = v0
      a = afield(x)
      t = 0
      while(t<tmax) {
        ddt = min(dt,tmax-t)
        v = v+a*ddt/2
        x = x+v*ddt
        a = afield(x)
        v = v+a*ddt/2
        t = t+ddt
      }
      return(list(x=x,v=v))
    }
    
    if (ti<=t0) {
      
      # leapfrog extrapolation before t0
      out = leapfrog(x0,-v0,t0-ti)
      x = out$x; v = out$v
      
    } else if (ti>=t1) {
      
      # leapfrog extrapolation after t1
      out = leapfrog(x1,v1,ti-t1)
      x = out$x; v = out$v
      
    } else {
      
      # leapfrog interpolation between t0 and t1
      out0 = leapfrog(x0,v0,ti-t0)
      out1 = leapfrog(x1,-v1,t1-ti)
      f = (t1-ti)/(t1-t0)
      x = f*out0$x+(1-f)*out1$x
      v = f*out0$v+(1-f)*out1$v
      
    }
    
  }
    
  return(list(x=x, v=v))
    
}


