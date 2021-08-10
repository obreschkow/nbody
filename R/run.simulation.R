#' @title Run a direct N-body simulation
#'
#' @description Run direct N-body simulations using a common adaptive time step.
#'
#' @importFrom Rcpp cppFunction
#' @importFrom utils read.table write.table
#'
#' @param sim structured list of simulation settings, which must contain the following sublists:\cr\cr
#'
#' \code{ics} is the sublist of initial conditions. It must contain the items. By default, they are all assumed in SI units, but if you can use any other units, as long as you specify the corresponding value of the gravitational constant \code{G} in the \code{para} sublist (see below).\cr
#' \code{m} = N-vector with the masses of the N particles\cr
#' \code{x} = N-by-3 matrix specifying the initial position in Cartesian coordinates\cr
#' \code{v} = N-by-3 matrix specifying the initial velocities\cr\cr
#'
#' \code{para} is the sublist of optional simulation parameters. It contains the items:\cr
#' \code{t.max} = final simulation time. If not given, a characteristic time is computed as \code{t.max = 2*pi*sqrt(R^3/GM)}, where \code{R} is the RMS radius and \code{M} is the total mass.\cr
#' \code{dt.max} = maximum time step. If not given, no maximum time step is imposed, meaning that the maximum time step is either equal to \code{dt.out} or the adaptive time step, whichever is smaller.\cr
#' \code{dt.min} = minimum time step used, unless a smaller time step is required to save an output or to land precisely on the final time \code{t.max}.
#' \code{dt.out} = output time step, i.e. time step between successive snapshots in the \code{output} sublist returned by \code{run.simulation}. If not given, \code{dt.max=t.max/100} is assumed.\cr
#' \code{eta} = scaling of adaptive time step. Smaller values lead to proportionally smaller adaptive time steps. Typical values range between 0.001 and 0.1. If not given, a default value of 0.01 is assumed. To use fixed time steps, set \code{eta=1e99} and set a time step \code{dt.max}.\cr
#' \code{integrator} = character string specifying the integrator to be used. Currently implemented integrators are 'euler' (1st order), 'leapfrog' (2nd order), 'yoshida' (4th order), 'yoshida6' (6th order). If not given, 'leapfrog' is the default integrator.\cr
#' \code{rsmooth} = optional smoothing radius. If not given, no smoothing is assumed.\cr
#' \code{afield} = a function(x,t) of positions \code{x} (N-by-3 matrix) and time \code{t} (scalar), specifying the external acceleration field. It must return an N-by-3 matrix. If not given, no external field is assumed.\cr
#' \code{G} = gravitational constant. If not given, the SI value specified in \code{\link{cst}} is used.
#'
#' \code{code} is the sublist of optional variables used when calling external simulation codes. It contains the items:\cr
#' \code{name} = character string specifying the name of the code, currently available options are "R" (default) and "nbodyx" (a simple, but fast N-body simulator in Fortran).
#' \code{file} = character string specifying the path+filename of the external compiled simulation code.
#' \code{interface} = character string specifying a temporary working path used as interface with external codes. If not given, the current working directory is used by default.
#'
#' @param measure.time logical flag that determines whether time computation time will be measured and displayed.
#'
#' @details
#' For a detailed description of the simulation method, please refer to the lecture notes on N-body simulations by Obreschkow (2019).
#'
#' @return The routine returns the structured list of the input argument, with one sublist \code{output} added. This sublist contains the items:
#' \item{t}{k-vector with the simulation times of the k snapshots.}
#' \item{x}{k-by-N-by-3 array giving the 3D coordinates of the N particles in k snapshots.}
#' \item{v}{k-by-N-by-3 array giving the 3D velocities of the N particles in k snapshots.}
#' \item{n.snapshots}{total number of snapshots.}
#' \item{n.iterations}{total number of iterations used to run the simulation.}
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
    message('\n')
    sim$output = NA
  }, finally = {
    return(sim)
  }))
}

# Core simulation routine
.run.sim = function(sim=NULL) {

  # check ICs
  if (is.null(sim)) stop('Mandatory argument "sim" is missing.\n')
  if (is.null(sim$ics)) stop('The initial conditions sublist "ics" is missing in the input argument sim.\n')
  if (length(sim$ics)!=3) stop('The sublist ics must contain exactly three items.\n')
  if (is.null(sim$ics$m)) stop('The item "ics$m" is missing.\n')
  if (is.null(sim$ics$x)) stop('The item "ics$x" is missing.\n')
  if (is.null(sim$ics$v)) stop('The item "ics$v" is missing.\n')
  sim$ics$x = as.matrix(sim$ics$x)
  sim$ics$v = as.matrix(sim$ics$v)
  if (length(sim$ics$x)!=3*length(sim$ics$m)) stop('Number of elements in x inconsistent with number of elements in m.\n')
  if (length(sim$ics$v)!=3*length(sim$ics$m)) stop('Number of elements in v inconsistent with number of elements in m.\n')
  if (length(sim$ics$m)==1) {
    sim$ics$x = rbind(as.vector(sim$ics$x))
    sim$ics$v = rbind(as.vector(sim$ics$v))
  }
  if (dim(sim$ics$x)[2]!=3) stop('x must be a matrix with 3 columns.\n')
  if (dim(sim$ics$v)[2]!=3) stop('v must be a matrix with 3 columns.\n')
  if (dim(sim$ics$x)[1]!=dim(sim$ics$v)[1]) stop('x and v must have the same number of rows.\n')
  if (dim(sim$ics$x)[1]!=length(sim$ics$m)) stop('The length of m must be equal to the number of rows of x.\n')
  if (dim(sim$ics$v)[1]!=length(sim$ics$m)) stop('The length of m must be equal to the number of rows of v.\n')

  # handle optional parameters
  if (is.null(sim$para)) sim$para = list()
  if (is.null(sim$para$G)) sim$para$G = cst$G
  if (is.null(sim$para$t.max)) {
    if (sim$para$G==0) {
      sim$para$t.max = 1
    } else {
      M = sum(abs(sim$ics$m))
      x0 = colSums(sim$ics$x*sim$ics$m)/M # center of mass
      R = sqrt(mean((sim$ics$x[,1]-x0[1])^2)+mean((sim$ics$x[,2]-x0[2])^2)+mean((sim$ics$x[,3]-x0[3])^2)) # RMS radius
      sim$para$t.max = 2*pi*sqrt(R^3/M/sim$para$G) # dynamical time scale
    }
  }
  if (is.null(sim$para$dt.max)) sim$para$dt.max=sim$para$t.max
  if (is.null(sim$para$dt.min)) sim$para$dt.min=0
  if (is.null(sim$para$dt.out)) sim$para$dt.out=sim$para$t.max/100
  if (is.null(sim$para$eta)) sim$para$eta=0.01
  if (is.null(sim$para$integrator)) sim$para$integrator='leapfrog'
  if (is.null(sim$para$rsmooth)) sim$para$rsmooth=0
  if (sim$para$t.max/sim$para$dt.out>1e5) stop('dt.out is too small for the simulation time t.max.\n')
  if (!is.null(sim$para$afield)) {
    a = try(sim$para$afield(sim$ics$x,0),silent=TRUE)
    if (length(dim(a))!=2) stop('afield is not a correctly vectorized function of (x,t).\n')
    if (any(dim(a)!=dim(sim$ics$x))) stop('afield is not a correctly vectorized function of (x,t).\n')
  }

  # handle external code
  if (is.null(sim$code)) sim$code = list()
  if (is.null(sim$code$name)) sim$code$name = 'R'
  if (sim$code$name!='R') {
    if (!is.null(sim$para$afield)) stop('afield can only be specified with code "R"')
    if (is.null(sim$code$file)) stop('code$file must be specified for external simulation codes')
    if (!file.exists(sim$code$file)) stop(sprintf('simulation code does not exist: %s',sim$code$file))
    if (file.access(sim$code$file,1)!=0) stop(sprintf('you have no permission to execute the code: %s',sim$code$file))
    if (is.null(sim$code$interface)) {
      if (file.access(getwd(),4)!=0) stop('you do not have permission to read in the current directory; please specify code$interface.')
      if (file.access(getwd(),2)!=0) stop('you do not have permission to write in the current directory; please specify code$interface.')
      sim$code$interface=file.path(getwd(),'tmp_nbody_interface')
    } else {
      if (file.access(sim$code$interface,4)!=0) stop(sprintf('you do not have permission to read in: %s',sim$code$interface))
      if (file.access(sim$code$interface,2)!=0) stop(sprintf('you do not have permission to write in: %s',sim$code$interface))
    }
    system(sprintf('rm -r %s',sim$code$interface))
    system(sprintf('mkdir -p %s',sim$code$interface))
    if (file.access(sim$code$interface,6)!=0) stop(sprintf('unable to create directory: %s',sim$code$interface))
  }

  # make global variables from ICs
  m = sim$ics$m
  n = length(m)
  x = sim$ics$x
  v = sim$ics$v


  # nbodyx (Fortran code by D. Obreschkow) ************************************************************************************

  if (sim$code$name=='nbodyx') {

    # make file names+paths
    filename.ics = file.path(sim$code$interface,'ics.txt')
    filename.para = file.path(sim$code$interface,'parameters.txt')
    filename.stats = file.path(sim$code$interface,'statistics.txt')
    filename.output = file.path(sim$code$interface,'snapshot_all.bin')

    # write ICs file
    ics = cbind(m,x,v)
    write.table(ics, filename.ics, col.names = FALSE, row.names = FALSE)

    # write parameter file
    para = {}
    para[1] = sprintf('inputfile %s',filename.ics)
    para[2] = sprintf('outputpath %s',sim$code$interface)
    para[3] = 'outputformat 3'
    para[4] = 'tinitial 0'
    para[5] = sprintf('tfinal %.15e',sim$para$t.max)
    para[6] = sprintf('dtmax %.15e',sim$para$dt.max)
    para[7] = sprintf('dtmin %.15e',sim$para$dt.min)
    para[8] = sprintf('dtout %.15e',sim$para$dt.out)
    para[9] = sprintf('G %.15e',sim$para$G)
    para[10] = sprintf('smoothing_radius %.15e',sim$para$rsmooth)
    para[11] = sprintf('eta %.15e',sim$para$eta)
    para[12] = sprintf('integrator %s',sim$para$integrator)
    write.table(para, file=filename.para, col.names = FALSE, row.names = FALSE, quote = FALSE)

    # run simulation
    system(sprintf('%s -parameterfile %s',sim$code$file, filename.para),intern=F)

    # read stats file
    stats = read.table(filename.stats)
    n.snapshots = stats[4,2]
    n.iterations = stats[5,2]
    n.acceleration.evaluations = stats[6,2]

    # read snapshot files
    if (!file.exists(filename.output)) stop(paste0('file does not exist: ',filename.output))
    ncheck = readBin(filename.output,'int',1,8)
    if (ncheck!=n) stop('wrong number of particles in file')
    dat = readBin(filename.output,'numeric',1+n.snapshots*(6*n+1),8)
    t.out = rep(NA,n.snapshots)
    x.out = v.out = array(NA,c(n.snapshots,n,3))
    for (i.out in seq(n.snapshots)) {
      offset = (6*n+1)*(i.out-1)
      t.out[i.out] = dat[offset+2]
      x.out[i.out,,] = dat[(offset+3):(offset+2+3*n)]
      v.out[i.out,,] = dat[(offset+3+3*n):(offset+2+6*n)]
    }

    # complete output
    sim$output = list(t = t.out[1:i.out], x = x.out[1:i.out,,], v = v.out[1:i.out,,],
                      n.snapshots = n.snapshots, n.iterations = n.iterations,
                      n.acceleration.evaluations = n.acceleration.evaluations)

  } else if (sim$code$name=='nbody6') {

      stop('nbody6 interface not yet available')

      # make file names+paths
      filename.input = file.path(sim$code$interface,'input')
      filename.output = file.path(sim$code$interface,'output')
      filename.ics = file.path(sim$code$interface,'fort.10')

      # write ICs file
      if (any(m<0)) stop('masses cannot be fixed (negative) with nbody6')
      ics = cbind(m,x,v)
      write.table(ics, filename.ics, col.names = FALSE, row.names = FALSE)

      # make parameter file (see details in nbody6/Ncode/define.f and nbody6/user/input_file_structure.pdf)
      lines = {}

      KSTART = 1 # Control index (1: new run; >1: restart; 3, 4, 5: new params).
      TCOMP = 1e3 # Maximum CPU time in minutes (saved in CPU).
      lines[1] = sprintf('%d %e',KSTART,TCOMP)

      N = n # Number of objects (N_s + 2*N_b; singles + 3*NBIN0 < NMAX).
      NFIX = 0 # ? Output frequency of data save or binaries (options 3 & 6).
      NCRIT = 0 # ? Final particle number (alternative termination criterion).
      NRAND = 0 # Random number seed.
      NNBMAX = 1 # ? Maximum number of neighbours (< LMAX - 5).
      NRUN = 1 # ? Run identification index.
      lines[2] = sprintf('%d %d %d %d %d %d',N,NFIX,NCRIT,NRAND,NNBMAX,NRUN)

      ETAI = 0.02 # Time-step parameter for irregular force polynomial.
      ETAR = 0.02 # Time-step parameter for regular force polynomial.
      RS0 = 0.0 # Initial radius of neighbour sphere (N-body units).
      DTADJ = 0.5 # Time interval for parameter adjustment (N-body units).
      DELTAT = 0.02 # ? Output time interval (N-body units).
      TCRIT = 1.0 # ? Termination time (N-body units).
      QE = 1e-3 # Energy tolerance (restart if DE/E > 5*QE & KZ(2) > 1).
      RBAR = 0.0 # Virial cluster radius in pc (set = 1 for isolated cluster).
      ZMBAR = 0.0 # Mean mass in solar units (=1.0 if 0; final depends on #20).
      lines[3] = sprintf('%e %e %e %e %e %e %e %e %e',ETAI,ETAR,RS0,DTADJ,DELTAT,TCRIT,QE,RBAR,ZMBAR)

      OPTIONS = rep(0,50)
      OPTIONS[1] = 0 # COMMON save unit 1 (=1: 'touch STOP'; =2: every 100*NMAX steps).
      OPTIONS[2] = 0 # COMMON save unit 2 (=1: at output; =2: restart if DE/E > 5*QE).
      OPTIONS[3] = 1 # Basic data unit 3 at output time
      OPTIONS[6] = 3
      OPTIONS[15] = 0 # ? Triple, quad, chain (#30 > 0) or merger search (>1: more output).
      OPTIONS[22] = 3 # ICs type: =3: no scaling of input read on fort.10;
      OPTIONS[26] = 0 # ? Slow-down of two-body motion (>=1: KS; >=2: chain; =3: rectify).
      OPTIONS[30] = 0 # Multiple regularization (=1: all; >1: BEGIN/END; >2: each step); =-1: CHAIN only; =-2: TRIPLE & QUAD only.
      OPTIONS[38] = 0 # Force polynomial corrections (=0: standard, no corrections; =1: all gains & losses included; =2: small FREG change skipped; =3: fast neighbour loss only).
      OPTIONS[40] = 0 # Neighbour number control (=1: increase if <NNB>  <  NNBMAX/2); >=2: fine-tuning at NNBMAX/5; =3: reduction of NNBMAX.

      lines[4] = gsub(',','',toString(OPTIONS[01:10]))
      lines[5] = gsub(',','',toString(OPTIONS[11:20]))
      lines[6] = gsub(',','',toString(OPTIONS[21:30]))
      lines[7] = gsub(',','',toString(OPTIONS[31:40]))
      lines[8] = gsub(',','',toString(OPTIONS[41:50]))

      DTMIN = 4e-5 # Time-step criterion for regularization search.
      RMIN = 0.0 # Distance criterion for regularization search.
      ETAU = 0.002 # Regularized time-step parameter (6.28/ETAU steps/orbit).
      ECLOSE = 0.0 # Binding energy per unit mass for hard binary (positive).
      GMIN = 1e-6 # Relative two-body perturbation for unperturbed motion.
      GMAX = 1e-3 # Secondary termination parameter for soft KS binaries.
      lines[9] = sprintf('%e %e %e %e %e %e',DTMIN,RMIN,ETAU,ECLOSE,GMIN,GMAX)

      ALPHAS = 0.0 # Power-law index for initial mass function (used if #20 < 2).
      BODY1 = 0.0 # Maximum particle mass before scaling (KZ(20): solar mass).
      BODYN = 0.0 # Minimum particle mass before scaling.
      NBIN0 = 0 # Number of primordial binaries (for IMF2 with KZ(20) > 1).
      NHI0 = 0 # Primordial hierarchies (may be needed in IMF if > 0).
      ZMET = 0.0 # Metal abundance (in range 0.03 - 0.0001).
      EPOCH0 = 0 # Evolutionary epoch (in 10**6 yrs; NB! < 0 for PM evolution).
      DTPLOT = 0 # Plotting interval for HRDIAG (N-body units; >= DELTAT).
      lines[10] = sprintf('%e %e %e %d %d %e %e %e',ALPHAS,BODY1,BODYN,NBIN0,NHI0,ZMET,EPOCH0,DTPLOT)

      Q = 0.0 # Virial ratio (Q = 0.5 for equilibrium).
      VXROT = 0 # XY-velocity scaling factor (> 0 for solid-body rotation).
      VZROT = 0 # Z-velocity scaling factor (not used if VXROT = 0).
      RTIDE = 0 # Unscaled tidal radius (#14 >= 2; otherwise copied to RSPH2).
      SMAX = 1 # Maximum time-step (factor of 2 commensurate with 1.0).
      lines[11] = sprintf('%e %e %e %e %e',Q,VXROT,VZROT,RTIDE,SMAX)
      write.table(lines, file=filename.input, col.names = FALSE, row.names = FALSE, quote = FALSE)

      # run simulation
      print(sprintf('cd %s; %s <input> output',sim$code$interface, sim$code$file))
      stop()
      system(sprintf('cd %s; %s <input> output',sim$code$interface, sim$code$file),intern=F)

  } else if (sim$code$name=='R') {

    if (sim$para$integrator=='euler') {

      iteration = function(dt) {
        .evaluate.accelerations()
        x <<- x+v*dt+0.5*a*dt^2
        v <<- v+a*dt
      }

    } else if (sim$para$integrator=='leapfrog') {

      iteration = function(dt) {
        v <<- v+a*dt/2
        x <<- x+v*dt
        .evaluate.accelerations()
        v <<- v+a*dt/2
      }

    } else if (sim$para$integrator%in%c('yoshida','yoshida4')) {

      .w1 = 1/(2-2^(1/3))
      .w0 = -2^(1/3)*.w1
      .c.yoshida = c(.w1/2,(.w0+.w1)/2,(.w0+.w1)/2,.w1/2)
      .d.yoshida = c(.w1,.w0,.w1)
      iteration = function(dt) {
        for (i in seq(3)) {
          x <<- x+.c.yoshida[i]*v*dt
          .evaluate.accelerations()
          v <<- v+.d.yoshida[i]*a*dt
        }
        x <<- x+.c.yoshida[4]*v*dt
      }

    } else if (sim$para$integrator=='yoshida6') {

      .x1 = 1/(2-2^(1/3))
      .x0 = -2^(1/3)*.x1
      .y1 = 1/(2-2^(1/5))
      .y0 = -2^(1/5)*.y1
      .d.yoshida = c(.x1*.y1,.x0*.y1,.x1*.y1,.x1*.y0,.x0*.y0,.x1*.y0,.x1*.y1,.x0*.y1,.x1*.y1)
      .c.yoshida = c(.d.yoshida[1]/2,(.d.yoshida[1:8]+.d.yoshida[2:9])/2,.d.yoshida[9]/2)
      iteration = function(dt) {
        for (i in seq(9)) {
          x <<- x+.c.yoshida[i]*v*dt
          .evaluate.accelerations()
          v <<- v+.d.yoshida[i]*a*dt
        }
        x <<- x+.c.yoshida[10]*v*dt
      }

    } else {

      stop('unknown integrator for this code')

    }

    # acceleration function
    .evaluate.accelerations = function() {
      n.acceleration.evaluations <<- n.acceleration.evaluations+1
      if (is.null(sim$para$afield)) {
        a[,] <<- 0
      } else {
        a = sim$para$afield(x,t)
      }
      f = accelerations(m,x,v,a,sim$para$G,rsmoothsqr)
      a[,] <<- f$a
      dt.var <<- f$dtvar
    }

    # initialize output variables
    n.out = ceiling(sim$para$t.max/sim$para$dt.out)+2
    t.out = rep(NA,n.out)
    x.out = v.out = array(NA,c(n.out,n,3))
    i.out = 0
    .save.snapshot = function() {
      i.out <<- i.out+1
      t.out[i.out] <<- t
      x.out[i.out,,] <<- x
      v.out[i.out,,] <<- v
    }

    # initialize first iteration
    dt.var = NULL # only used for variable time-stepping
    t = 0
    n.iterations = 0
    n.acceleration.evaluations = 0
    .save.snapshot()
    t.next = sim$para$dt.out # time of next output
    rsmoothsqr = sim$para$rsmooth^2
    a = array(0,c(n,3)) # this line is necessary to make 'a' a local variable
    .evaluate.accelerations()

    # iterate
    while (t < sim$para$t.max) {
      dt = min(sim$para$dt.max, sim$para$t.max-t, t.next-t, max(sim$para$dt.min, sim$para$eta*dt.var))
      iteration(dt)
      t = t+dt
      n.iterations = n.iterations+1
      if (t>=t.next | t>=sim$para$t.max) {
        .save.snapshot()
        t.next = i.out*sim$para$dt.out
      }
    }

    # finalise output
    sim$output = list(t = t.out[1:i.out], x = x.out[1:i.out,,], v = v.out[1:i.out,,],
                      n.snapshots = i.out, n.iterations = n.iterations, n.acceleration.evaluations = n.acceleration.evaluations)

  } else {

    stop(sprintf('unknown simulation code: %s',sim$code$name))

  }

  # turn list into a simulation class
  class(sim) = 'simulation'

  return(sim)

}
