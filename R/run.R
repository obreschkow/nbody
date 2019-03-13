#' Run a direct N-body simulation
#'
#' This function finds the most likely P-dimensional model parameters of a D-dimensional distribution function (DF) generating an observed set of N objects with D-dimensional observables x, accounting for measurement uncertainties and a user-defined selection function. For instance, if the objects are galaxies, \code{dffit} can fit a mass function (D=1), a mass-size distribution (D=2) or the mass-spin-morphology distribution (D=3). A full description of the algorithm can be found in Obreschkow et al. (2017).
#'
#' @importFrom Rcpp cppFunction
#'
#' @param sim is a list of simulation settings:\cr\cr
#' \code{x} is a M-by-D array of M points, defining a regular cartesian grid in the D-dimensional observable space.\cr\cr
#' \code{xmin} and \code{xmax} are D-vectors specifying the lower and upper boundary of the grid in the D-dimensional observable space.\cr\cr
#' \code{dx} is a D-vector specifying the steps between grid points.\cr\cr
#'
#' @param measure.time is a logical flag that determines whether time computation time will be measured and displayed.
#'
#' @details
#' For a detailed description of the method, please refer to the peer-reviewed publication by Obreschkow et al. 2017 (in prep.).
#'
#' @return The routine returns a structured list \code{sim} with the following sublists:
#'
#' \item{selection}{is a list describing the selection function of the data used when fitting the generative DF. The most important entries in this list are:\cr\cr
#' \code{veff} is a function of a D-dimensional vector, specifying the effective volume associated with an object of properties xval. If LSS is corrected for, i.e. if the argument \code{correct.lss.bias} is set to \code{TRUE}, this function is the final effecive volume, including the effect of LSS.\cr\cr
#' \code{veff.no.lss} is a function of a D-dimensional vector, specifying the effecive volume associate with an object, if LSS were not accounted for. This function is indentical to \code{veff}, if LSS is not accounted for in the fit, i.e. if the argument \code{correct.lss.bias} is set to \code{FALSE}.\cr\cr
#' \code{mode} is an integer specifying the format of the input argument \code{selection}.}
#'
#' \item{fit}{is a list describing the fitted generative distribution function. Its most important entries are:\cr\cr
#' \code{p.best} is a P-vector giving the most likely model parameters according to the MML method.\cr\cr
#' \code{p.sigma} and \code{p.covariance} are the standard deviations and covariance matrices of the best-fitting parameters in the Gaussian approximation from the Hessian matrix of the modified likelihood function.\cr\cr
#'  \code{gdf} is a function of a D-dimensional vector, which is the generative DF, evaluated at the parameters \code{p.best}.\cr\cr
#'  \code{scd} is a function of a D-dimensional vector, which gives the predicted source counts of the most likely model, i.e. scd(x)=gdf(x)*veff(x).}
#'
#' @keywords N-body simulation
#'
#' @examples
#' sim = list(ics = ics.halley(), para = list(t.max = 76*cst$yr, dt.out = 3*cst$month))
#' sim = run.simulation(sim)
#' plot.simulation(sim, cst$AU, xlim=c(-20,60), ylim=c(-40,40), xlab='[AU]', ylab='[AU]')
#'
#' @author Danail Obreschkow
#'
#' @export run.simulation

####################################################################################
# simulation code ##################################################################
####################################################################################

# Wrapper to catch errors and measure simulation time
run.simulation = function(sim, measure.time = TRUE) {
  return(tryCatch({
    t.start = Sys.time()
    sim = .run.sim(sim)
    if (measure.time) {
      t.end = Sys.time()
      time.taken = as.double(t.end)-as.double(t.start)
      message(sprintf('Simulation successfully completed in %.2fs.',time.taken))
    }
  }, warning = function(w) {
    message('Warning(s) produced in executing the simulation.')
    sim$output = NA
  }, error = function(e) {
    message('Error(s) produced in executing the simulation.')
    message(e)
    sim$output = NA
  }, finally = {
    return(sim)
  }))
}

# Core simulation routine
.run.sim = function(sim) {

  # integration algorithms
  .iteration = {}
  .iteration$euler = function(dt) {
    .evaluate.accelerations()
    x <<- x+v*dt+0.5*a*dt^2
    v <<- v+a*dt
  }

  .iteration$leapfrog = function(dt) {
    v <<- v+a*dt/2
    x <<- x+v*dt
    .evaluate.accelerations()
    v <<- v+a*dt/2
  }

  .w1 = 1/(2-2^(1/3))
  .w0 = -2^(1/3)*.w1
  .c.yoshida = c(.w1/2,(.w0+.w1)/2,(.w0+.w1)/2,.w1/2)
  .d.yoshida = c(.w1,.w0,.w1)
  .iteration$yoshida = function(dt) {
    for (i in seq(3)) {
      x <<- x+.c.yoshida[i]*v*dt
      .evaluate.accelerations()
      v <<- v+.d.yoshida[i]*a*dt
    }
    x <<- x+.c.yoshida[4]*v*dt
  }

  # acceleration function
  .evaluate.accelerations = function() {
    a[,] <<- 0
    if (length(rsmoothsqr)==0) rsmoothsqr = 0
    f = accelerations(m,x,cst$G,rsmoothsqr)
    a[,] <<- f$a
    dt.var <<- f$dtvar
  }

  # make global variables from ICs
  m = sim$ics$m
  x = sim$ics$x
  v = sim$ics$v
  a = array(0,c(3,n))
  n = length(m)

  # set iteration function
  if (is.null(sim$para$algorithm)) sim$para$algorithm = 'leapfrog'
  custom.iteration = .iteration[[sim$para$algorithm]]

  # initialize output variables
  n.out = floor(sim$para$t.max/sim$para$dt.out)+2
  x.out = v.out = array(NA,c(n.out,3,n))
  i.out = 1
  x.out[1,,] = x
  v.out[1,,] = v
  t.out = sim$para$dt.out # time of next output

  # advanced settings (for future use)
  if (is.null(sim$para$rsmooth)) {
    rsmoothsqr = 0
  } else {
    rsmoothsqr = sim$para$rsmooth^2 # only used for smoothed particles
  }
  if (is.null(sim$para$eta)) sim$para$eta = 0.01
  dt.var = NULL # only used for variable time-stepping

  # prepare first iteration
  t = 0
  n.iterations = 0
  .evaluate.accelerations()

  # iterate
  while (t < sim$para$t.max) {
    dt = min(sim$para$dt.max, sim$para$t.max-t, t.out-t, sim$para$eta*dt.var)
    custom.iteration(dt)
    t = t+dt
    if (t>=t.out) {
      i.out = i.out+1
      x.out[i.out,,] = x
      v.out[i.out,,] = v
      t.out = t.out+sim$para$dt.out
    }
    n.iterations = n.iterations+1
  }

  # finalise output
  i.out = i.out+1
  x.out[i.out,,] = x
  v.out[i.out,,] = v
  sim$output = list(x.out = x.out[1:i.out,,], v.out = v.out[1:i.out,,],
                    n.iterations = n.iterations)
  return(sim)

}
