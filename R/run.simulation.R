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
#' \code{para} is the sublist of optional simulation parameters. It contains the items:\cr
#' \code{t.max} = final simulation time. If not given, a characteristic time is computed as \code{t.max = 2*pi*sqrt(R^3/GM)}, where \code{R} is the RMS radius and \code{M} is the total mass.\cr
#' \code{dt.max} = maximum time step. If not given, no maximum time step is imposed, meaning that the maximum time step is either equal to \code{dt.out} or the adaptive time step, whichever is smaller.\cr
#' \code{dt.min} = minimum time step used, unless a smaller time step is required to save an output or to land precisely on the final time \code{t.max}.
#' \code{dt.out} = output time step, i.e. time step between successive snapshots in the \code{output} sublist returned by \code{run.simulation}. If not given, \code{dt.max=t.max/100} is assumed.\cr
#' \code{eta} = scaling of adaptive time step. Smaller values lead to proportionally smaller adaptive time steps. Typical values range between 0.001 and 0.1. If not given, a default value of 0.01 is assumed. To use fixed time steps, set \code{eta=1e99} and set a time step \code{dt.max}.\cr
#' \code{algorithm} = character string specifying the integrator to be used. Currently implemented algorithms are 'euler' (1st order), 'leapfrog' (2nd order), 'yoshida' (4th order). If not given, 'leapfrog' is the default algorithm.\cr
#' \code{rsmooth} = optional smoothing radius. If not given, no smoothing is assumed.\cr
#' \code{afield} = a function(x,t) of position \code{x} (3-vector) and time \code{t} (scalar), specifying the external acceleration field. If not given, no external field is assumed.
#' \code{G} = alternative gravitational constant. If not given, the true value specified in \code{\link{cst}} is used.
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
.run.sim = function(sim=NULL) {

  # check ICs
  if (is.null(sim)) stop('Mandatory argument "sim" is missing.')
  if (is.null(sim$ics)) stop('The initial conditions sublist "ics" is missing in the input argument sim.')
  if (length(sim$ics)!=3) stop('The sublist ics must contain exactly three items.')
  if (is.null(sim$ics$m)) stop('The item "ics$m" is missing.')
  if (is.null(sim$ics$x)) stop('The item "ics$x" is missing.')
  if (is.null(sim$ics$v)) stop('The item "ics$v" is missing.')
  if (dim(sim$ics$x)[2]!=3) stop('x must be a matrix with 3 columns.')
  if (dim(sim$ics$v)[2]!=3) stop('v must be a matrix with 3 columns.')
  if (dim(sim$ics$x)[1]!=dim(sim$ics$v)[1]) stop('x and v must have the same number of rows.')
  if (dim(sim$ics$x)[1]!=length(sim$ics$m)) stop('The length of m must be equal to the number of rows of x.')
  if (dim(sim$ics$v)[1]!=length(sim$ics$m)) stop('The length of m must be equal to the number of rows of v.')

  # handle optional parameters
  if (is.null(sim$para)) sim$para = list()
  if (is.null(sim$para$G)) sim$para$G = cst$G
  if (is.null(sim$para$t.max)) {
    M = sum(sim$ics$m)
    x0 = colSums(sim$ics$x*sim$ics$m)/M # center of mass
    R = sqrt(mean((sim$ics$x[,1]-x0[1])^2)+mean((sim$ics$x[,2]-x0[2])^2)+mean((sim$ics$x[,3]-x0[3])^2)) # RMS radius
    sim$para$t.max = 2*pi*sqrt(R^3/M/sim$para$G) # dynamical time scale
  }
  if (is.null(sim$para$dt.out)) sim$para$dt.out=sim$para$t.max/100
  if (is.null(sim$para$eta)) sim$para$eta=0.01
  if (is.null(sim$para$algorithm)) sim$para$algorithm='leapfrog'
  if (is.null(sim$para$rsmooth)) sim$para$rsmooth=0
  if (sim$para$t.max/sim$para$dt.out>1e4) stop('dt.out is too small for the simulation time t.max.')

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

  # prepare first iteration
  dt.var = NULL # only used for variable time-stepping
  custom.iteration = .iteration[[sim$para$algorithm]]
  rsmoothsqr = sim$para$rsmooth^2
  t = 0
  n.iterations = 0
  .evaluate.accelerations()

  # iterate
  while (t < sim$para$t.max) {
    dt = min(sim$para$dt.max, sim$para$t.max-t, t.next-t, max(sim$para$dt.min, sim$para$eta*dt.var))
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
