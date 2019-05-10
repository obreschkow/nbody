#' @title Run a direct N-body simulation
#'
#' @description Run direct N-body simulations using a common adaptive time step.
#'
#' @importFrom Rcpp cppFunction
#'
#' @param sim structured list of simulation settings, which must contain the following sublists:\cr\cr
#'
#' \code{ics} is the sublist of initial conditions. It must contain the items: \cr
#' \code{m} = N-vector with the masses of the N particles\cr
#' \code{x} = N-by-3 matrix specifying the initial position in cartesian coordinates\cr
#' \code{v} = N-by-3 matrix specifying the initial velocities\cr\cr
#'
#' \code{para} is the sublist of simulation parameters. It contains the items:\cr
#' \code{t.max} = final simulation time.\cr
#' \code{dt.xmax} = maximum time step; ignored if set to \code{NULL} (default).\cr
#' \code{dt.out} = output time step, i.e. time step between successive snapshots in the \code{output} sublist returned by \code{run.simulation}.\cr
#' \code{eta} = scaling of adaptive time step. Smaller values lead to proportionally smaller adaptive time steps. Typical values range between 0.001 and 0.1. The default is 0.01.\cr
#' \code{algorithm} = character string specifying the integrator to be used. Currently implemented algorithms are 'euler' (1st order), 'leapfrog' (2nd order, default), 'yoshida' (4th order).\cr
#' \code{rsmooth} = optional smoothing radius. No smoothing is assumed by default.\cr
#' \code{afield} = a function(x,t) of position \code{x} (3-vector) and time \code{t} (scalar), specifying the external acceleration field. No field is assumes if set to \code{NULL} (default).
#' \code{G} = alternative gravitational constant. If set to \code{NULL}, the true value specified in \code{\link{cst}} will be used.
#'
#' @param measure.time logical flag that determines whether time computation time will be measured and displayed.
#'
#' @details
#' For a detailed description of the simulation method, please refer to the lecture notes on N-body simulations by Obreschkow (2019).
#'
#' @return The routine returns the structured list of the input argument, with one sublist \code{output} added. This sublist contains the items:
#' \item{t}{k-vector with the simulation times of the k output snapshots.}
#' \item{x}{k-by-N-by-3 array giving the 3D coordinates of the N particles in k snapshots.}
#' \item{v}{k-by-N-by-3 array giving the 3D velocities of the N particles in k snapshots.}
#' \item{n.iterations}{Number of iterations used to run the simulation.}
#'
#' @keywords N-body simulation
#'
#' @examples
#' sim = setup.halley()
#' sim = run.simulation(sim)
#' plot(sim, units=cst$AU, xlim=c(-20,60), ylim=c(-40,40), xlab='[AU]', ylab='[AU]')
#' cat(sprintf('This simulation was run with %d iterations.\n',sim$output$n.iterations))
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

  # check ICs
  if (is.null(sim$ics)) stop('The initial conditions sublist "ics" is missing in the input argument sim.')
  if (length(sim$ics)!=3) stop('The sublist ics must contain exactly three items.')
  if (is.null(sim$ics$m)) stop('The item "ics$m" is missing.')
  if (is.null(sim$ics$x)) stop('The item "ics$x" is missing.')
  if (is.null(sim$ics$v)) stop('The item "ics$v" is missing.')

  # the parameters
  if (is.null(sim$para)) stop('The initial conditions sublist "para" is missing in the input argument sim.')
  if (is.null(sim$para$t.max)) stop('The item "para$t.max" is missing.')

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
    if (is.null(sim$para$afield)) {
      a[,] <<- 0
    } else {
      a <<- sim$para$afield(x,t)
    }
    if (length(rsmoothsqr)==0) rsmoothsqr = 0
    f = accelerations(m,x,a,sim$para$G,rsmoothsqr)
    a[,] <<- f$a
    dt.var <<- f$dtvar
  }


  # initialize parameters
  if (is.null(sim$para$G)) sim$para$G = cst$G
  if (is.null(sim$para$algorithm)) sim$para$algorithm = 'leapfrog'
  custom.iteration = .iteration[[sim$para$algorithm]]
  if (is.null(sim$para$rsmooth)) sim$para$rsmooth=0
  rsmoothsqr = sim$para$rsmooth^2
  if (is.null(sim$para$eta)) sim$para$eta = 0.01

  # make global variables from ICs
  m = sim$ics$m
  n = length(m)
  x = sim$ics$x
  v = sim$ics$v
  a = array(0,c(n,3)) # this line is necessary to make 'a' a local variable

  # initialize output variables
  n.out = ceiling(sim$para$t.max/sim$para$dt.out)+2
  t.out = rep(NA,n.out)
  x.out = v.out = array(NA,c(n.out,n,3))
  i.out = 1
  x.out[1,,] = x
  v.out[1,,] = v
  t.next = sim$para$dt.out # time of next output

  dt.var = NULL # only used for variable time-stepping

  # prepare first iteration
  t = 0
  n.iterations = 0
  .evaluate.accelerations()

  # iterate
  while (t < sim$para$t.max) {
    dt = min(sim$para$dt.max, sim$para$t.max-t, t.next-t, sim$para$eta*dt.var)
    custom.iteration(dt)
    t = t+dt
    n.iterations = n.iterations+1
    if (t>=t.next & t<sim$para$t.max) {
      i.out = i.out+1
      t.out[i.out] = t
      x.out[i.out,,] = x
      v.out[i.out,,] = v
      t.next = t.next+sim$para$dt.out
    }
  }

  # finalise output
  i.out = i.out+1
  t.out[i.out] = t
  x.out[i.out,,] = x
  v.out[i.out,,] = v
  sim$output = list(t = t.out[1:i.out], x = x.out[1:i.out,,], v = v.out[1:i.out,,],
                    n.iterations = n.iterations)

  # turn list into a simulation class
  class(sim) = 'simulation'

  return(sim)

}
