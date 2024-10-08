#' @title Run a direct N-body simulation
#'
#' @description Run direct N-body simulations using an adaptive block timestep.
#'
#' @importFrom Rcpp cppFunction
#' @importFrom utils read.table write.table
#' @importFrom stats sd
#'
#' @param sim structured list of simulation settings, which must contain the following sublists:\cr\cr
#'
#' \code{ics} is the sublist of initial conditions. It must contain the items:\cr
#' \code{m} = N-vector with the masses of the N particles. Negative mass values are considered as positive masses belonging to a background field, which is not subject to any forces. Therefore particles with negative mass will have a normal effect on particles with positive masses, but they will not, themselves, be accelerated by any other particle.\cr
#' \code{x} = N-by-3 matrix specifying the initial position in Cartesian coordinates\cr
#' \code{v} = N-by-3 matrix specifying the initial velocities\cr\cr
#'
#' \code{para} is an optional sublist of optional simulation parameters. It contains the items:\cr
#' \code{t.max} = final simulation time in simulation units (see details). If not given, a characteristic time is computed as \code{t.max = 2*pi*sqrt(R^3/GM)}, where \code{R} is the RMS radius and \code{M} is the total mass.\cr
#' \code{dt.max} = maximum time step. If not given, no maximum time step is imposed, meaning that the maximum time step is either equal to \code{dt.out} or the adaptive time step, whichever is smaller.\cr
#' \code{dt.min} = minimum time step used, unless a smaller time step is required to save an output or to land precisely on the final time \code{t.max}.
#' \code{dt.out} = output time step, i.e. time step between successive snapshots in the \code{output} sublist returned by \code{run.simulation}. If not given, \code{dt.max=t.max/100} is assumed.\cr
#' \code{eta} = accuracy parameter of adaptive time step. Smaller values lead to proportionally smaller adaptive time steps. Typical values range between 0.001 and 0.1. If not given, a default value of 0.01 is assumed. To use fixed time steps, set \code{eta=1e99} and set a time step \code{dt.max}.\cr
#' \code{integrator} = character string specifying the integrator to be used. Currently implemented integrators are 'euler' (1st order), 'leapfrog' (2nd order), 'yoshida' (4th order), 'yoshida6' (6th order). If not given, 'leapfrog' is the default integrator.\cr
#' \code{rsmooth} = optional smoothing radius. If not given, no smoothing is assumed.\cr
#' \code{afield} = a function(x,t) of positions \code{x} (N-by-3 matrix) and time \code{t} (scalar), specifying the external acceleration field. It must return an N-by-3 matrix. If not given, no external field is assumed. If the external code "nbodyx" is used, then afield should be a vector of the parameters p1, p2, ... for the external acceleration field of "nbodyx".\cr
#' \code{G} = gravitational constant in simulation units (see details). If not given, the measured value in SI units is used.\cr
#' \code{box.size} = scalar>=0. If 0, open boundary conditions are adopted. If >0, the simulation is run in a cubic box of side length box.size with periodic boundary conditions. In this case, the cubic box is contained in the interval [0,box.size) in all three Cartesian coordinates, and all initial positions must be contained in this interval. For periodic boundary conditions, the force between any two particles is always calculated along their shortest separation, which may cross 0-3 boundaries. The exception is GADGET-4, which also evaluates the forces from the periodic repetitions.\cr
#' \code{include.bg} = logical argument. If FALSE (default), only foreground particles, i.e. particles with masses >=0, are contained in the output vectors \code{x} and \code{v}. If TRUE, all particles are included.\cr
#' \code{ignore.size} = logical argument. If FALSE (default), the number of snapshots is limited to 1e7 and an error is produced if this number is expected to be exceeded.\cr\cr
#'
#' \code{code} is an optional sublist to force the use of an external simulation code (see details). It contains the items:\cr
#' \code{name} = character string specifying the name of the code, currently available options are "R" (default), "nbodyx" (a simple, but fast N-body simulator in Fortran) and "gadget4" (a powerful N-body+SPH simulator, not very adequate for small direct N-body simulations).\cr
#' \code{file} = character string specifying the path+filename of the external compiled simulation code.\cr
#' \code{interface} = optional character string specifying a temporary working path used as interface with external codes. NOTE: All existing files in this directory are deleted! If not given, the current working directory is used by default.\cr
#' \code{kind} = optional number of bytes per floating-point number used in nbodyx output files (has no bearing on computation accuracy)\cr
#' \code{gadget.np} = number of processors used with GADGET-4 (defaults to 1, which is normally best for small direct N-body runs)\cr
#'
#' @param measure.time logical flag that determines whether time computation time will be measured and displayed.
#'
#' @param verbose logical flag indicating whether to show console outputs from external codes. Ignored when using the in-built simulator.
#'
#' @details
#' UNITS: The initial conditions (in the sublist \code{ics}) can be provided in any units. The units of mass, length and velocity then fix the other units.
#' For instance, [unit of time in seconds] = [unit of length in meters] / [unit of velocity in m/s]. E.g., if initial positions are given in units of 1AU=1.49598e11m and velocities in units of 1km/s, one unit of time is 1.49598e8s=4.74yrs.
#' Likewise, units of the gravitational constant \code{G} are given via [unit of G in m^3*kg^(-1)*s^(-2)] = [unit of length in meters] * [unit of velocity in m/s]^2 / [unit of mass in kg]. E.g., for length units of 1AU=1.49598e11m, velocity units of 1km/s=1e3m/s and mass units of 1Msun=1.98847e30kg, a unit of G is
#' 7.523272e-14 m^3*kg^(-1)*s^(-2). In these units the true value of G is about 887.154.\cr\cr
#'
#' NBODYX simulator:\cr
#' Can be downloaded from github via\cr
#' \code{git clone https://github.com/obreschkow/nbodyx}\cr
#' Details on installing, compiling and running the code are given in the README file.\cr
#' Note: To run very high-accuracy simulations, such as the Pythagorean three-body problem, you can use 128-bit floating-point numbers by compiling the code as\cr
#' \code{make kind=16}\cr\cr
#'
#' GADGET-4 simulator:\cr
#' This his a very powerful N-body+SPH simulator used primarily for large astrophysical simulations. GADGET-4 is not particularly suitable for small direct N-body problems, but it can nonetheless be used for such simulations for the sake of comparison, at least if not too much accuracy is needed and if a massively increased computational overhead is acceptable.
#' Please refer to https://wwwmpa.mpa-garching.mpg.de/gadget4 for details on how to download and compile the code. In order to use GADGET-4 with this R-package, it must be compiled with the following compile-time options (in the file Config.sh):\cr
#' \code{NTYPES=2}\cr
#' \code{GADGET2_HEADER}\cr
#' \code{SELFGRAVITY}\cr
#' \code{ALLOW_DIRECT_SUMMATION}\cr
#' \code{HIERARCHICAL_GRAVITY}\cr
#' \code{DOUBLEPRECISION=1}\cr
#' \code{ENLARGE_DYNAMIC_RANGE_IN_TIME}\cr
#' If and only if periodic boundary conditions are used, you also need to add the option\cr
#' \code{PERIODIC}\cr
#' If you plan to often switch between runs with open and periodic boundaries, it may be advisable to compile two versions of GADGET-4, with and without this option. To do so, one needs to create two sub-directories with the respective Config.sh files and compile them via\cr
#' \code{make -j [number of cores] DIR=[path containing Config.sh with PERIODIC]}\cr
#' \code{make -j [number of cores] DIR=[path containing Config.sh without PERIODIC]}\cr
#' The runtime parameter file (param.txt) needed by GADGET-4 is written automatically when calling \code{run.simulation}. The gravitational softening length in GADGET-4 is computed as sim$para$rsmooth/2.8, which ensures that the particles behave like point masses at separations beyond sim$para$rsmooth. If rsmooth is not provided, it is computed as \code{stats::sd(apply(sim$ics$x,2,sd))*1e-5}. The accuracy parameter ErrTolIntAccuracy is set equal to sim$para$eta/sim$para$rsmooth*1e-3, which gives roughly comparable accuracy to in-built simulator for the Leapfrog integrator.
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
#' AU = 149597870700 # Astronomical unit in meters
#' plot(sim, units=AU, xlim=c(-20,60), ylim=c(-40,40), xlab='[AU]', ylab='[AU]')
#' cat(sprintf('This simulation was run with %d iterations.\n',sim$output$n.iterations))
#'
#' @author Danail Obreschkow
#'
#' @export run.simulation

####################################################################################
# simulation code ##################################################################
####################################################################################

# Wrapper to catch errors and measure simulation time
run.simulation = function(sim, measure.time=TRUE, verbose=TRUE) {
  tryCatch({
    t.start = Sys.time()
    sim = .run.sim(sim,verbose)
    if (measure.time) {
      t.end = Sys.time()
      time.taken = as.double(t.end)-as.double(t.start)
      message(sprintf('Simulation successfully completed in %.2fs.',time.taken))
    }
  }, warning = function(w) {
    message('Warning(s) produced in executing the simulation.')
    sim$output = NA
  }, error = function(e) {
    stop(paste('run.simulation: ',e[[1]]),call.=FALSE)
  })
  return(sim)
}

# Core simulation routine
.run.sim = function(sim=NULL,verbose) {

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
  if (is.null(sim$para$G)) sim$para$G = 6.67430e-11 # gravitational constant in SI units
  if (is.null(sim$para$box.size)) sim$para$box.size = 0
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
  if (sim$para$rsmooth<0) stop('smoothing radius cannot be negative.')
  if (sim$para$t.max/sim$para$dt.out>1e7) {
    if (is.null(sim$para$ignore.size)||!sim$para$ignore.size) {
      stop('dt.out is too small for the simulation time t.max\nSupress error by setting parameter ignore.size to true.\n')
    }
  }
  if (is.null(sim$para$include.bg)) sim$para$include.bg = FALSE

  # handle external code
  if (is.null(sim$code)) {
    if (is.null(.nbody.env$code)) {
      sim$code = list(name='R')
    } else {
      sim$code = .nbody.env$code
    }
  } else {
    if (is.null(sim$code$name)) stop('A code$name must be specified, if the argument code is given.')
  }
  if (sim$code$name!='R') {
    if (is.null(sim$code$file)) stop('code$file must be specified for external simulation codes')
    if (!file.exists(sim$code$file)) stop(sprintf('simulation code does not exist: %s',sim$code$file))
    if (file.access(sim$code$file,1)!=0) stop(sprintf('you have no permission to execute the code: %s',sim$code$file))
    if (is.null(sim$code$interface)) {
      sim$code$interface=tempfile()
      dir.create(sim$code$interface)
    }
    if (file.access(sim$code$interface,4)!=0) stop(sprintf('you do not have permission to read in: %s',sim$code$interface))
    if (file.access(sim$code$interface,2)!=0) stop(sprintf('you do not have permission to write in: %s',sim$code$interface))
    files = list.files(sim$code$interface, full.names = TRUE, recursive = TRUE)
    unlink(files,recursive=TRUE) # remove all files from interface directory
  }

  # handle external acceleration field
  if (!is.null(sim$para$afield)) {
    if (sim$code$name=='R') {
      a = try(sim$para$afield(sim$ics$x,0),silent=TRUE)
      if (length(dim(a))!=2) stop('afield is not a correctly vectorized function of (x,t).\n')
      if (any(dim(a)!=dim(sim$ics$x))) stop('afield is not a correctly vectorized function of (x,t).\n')
    } else if (sim$code$name=='nbodyx') {
      if (!is.numeric(sim$para$afield)) stop('When using nbodyx with an external field, then afield should be the vector of parameters. Set afield=0, if no parameters used.')
    } else {
      stop('afield can only be specified for codes "R" and "nbodyx"')
    }
  }

  # nbodyx (Fortran code by D. Obreschkow) ************************************************************************************

  if (sim$code$name=='nbodyx') {

    # handle input
    if (is.null(sim$code$kind)) sim$code$kind=8

    # make file names+paths
    filename.ics = file.path(sim$code$interface,'ics.txt')
    filename.para = file.path(sim$code$interface,'parameters.txt')
    filename.stats = file.path(sim$code$interface,'statistics.txt')
    filename.output = file.path(sim$code$interface,'snapshot_all.bin')

    # write ICs file
    ics = cbind(sim$ics$m,sim$ics$x,sim$ics$v)
    write.table(ics, filename.ics, col.names = FALSE, row.names = FALSE)

    # write parameter file
    para = {}
    para[1] = sprintf('inputfile %s',filename.ics)
    para[2] = sprintf('outputpath %s',sim$code$interface)
    para[3] = 'outputformat 3'
    para[4] = sprintf('kind %d',sim$code$kind)
    para[5] = sprintf('include_bg %d',as.numeric(sim$para$include.bg))
    para[6] = 'tinitial 0'
    para[7] = sprintf('tfinal %.15e',sim$para$t.max)
    para[8] = sprintf('dtmax %.15e',sim$para$dt.max)
    para[9] = sprintf('dtmin %.15e',sim$para$dt.min)
    para[10] = sprintf('dtout %.15e',sim$para$dt.out)
    para[11] = sprintf('G %.15e',sim$para$G)
    para[12] = sprintf('smoothing_radius %.15e',sim$para$rsmooth)
    para[13] = sprintf('eta %.15e',sim$para$eta)
    para[14] = sprintf('box_size %.15e',sim$para$box.size)
    para[15] = sprintf('integrator %s',sim$para$integrator)
    if (is.null(sim$para$afield)) {
      para[16] = sprintf('acceleration n')
    } else {
      para[16] = sprintf('acceleration y')
      for (i in seq_along(sim$para$afield)) {
        para[16+i] = sprintf('p%d %.15e',i,sim$para$afield[i])
      }
    }

    write.table(para, file=filename.para, col.names = FALSE, row.names = FALSE, quote = FALSE)

    # run simulation
    error.code = system(sprintf('%s -parameterfile %s',sim$code$file, filename.para),intern=FALSE,ignore.stderr=TRUE,ignore.stdout=!verbose)
    if (error.code>0) stop('Could not execute the nbodyx code. Try recompiling nbodyx.')

    # read stats file
    if (!file.exists(filename.stats)) stop(paste0('file does not exist: ',filename.stats))
    stats = read.table(filename.stats)
    n.snapshots = stats[4,2]
    n.iterations = stats[5,2]
    n.acceleration.evaluations = stats[6,2]

    # read output of nbodyx
    if (!file.exists(filename.output)) stop(paste0('file does not exist: ',filename.output))
    n.save = ifelse(sim$para$include.bg,length(sim$ics$m),sum(sim$ics$m>=0))
    if (n.save==0) stop('No particles in simulation output. Consider setting include.bg=TRUE.')
    fileid = file(filename.output, "rb")
    ncheck = readBin(fileid,'int',1,sim$code$kind)
    if (ncheck!=n.save) stop('wrong number of particles in file')
    dat = readBin(fileid,'numeric',1+n.snapshots*(6*n.save+1),sim$code$kind)
    close(fileid)

    # restructure output
    dat = array(dat,c(1+6*n.save,n.snapshots))
    t.out = dat[1,]
    x.out = drop(array(t(dat[2:(1+3*n.save),]),c(n.snapshots,n.save,3)))
    v.out = drop(array(t(dat[(2+3*n.save):(1+6*n.save),]),c(n.snapshots,n.save,3)))

    # complete output
    sim$output = list(t = t.out, x = x.out, v = v.out,
                      n.snapshots = n.snapshots, n.iterations = n.iterations,
                      n.acceleration.evaluations = n.acceleration.evaluations)

  }


  # GADGET-4 (Cosmological N-body simulator by V. Springel) ****************************************************************

  else if (sim$code$name=='gadget4') {

    # make file names+paths
    filename.ics = file.path(sim$code$interface,'ics.dat')
    filename.para = file.path(sim$code$interface,'param.txt')
    directory.output = file.path(sim$code$interface,'output')

    # write ICs file
    .gadget_write(sim$ics, filename.ics)

    # write parameter file
    if (sim$para$rsmooth==0) sim$para$rsmooth = stats::sd(apply(sim$ics$x,2,sd))*1e-5 # default smoothing radius if not provided
    smoothing = sim$para$rsmooth/2.8 # => the particles will behave exactly as point masses at separations larger than rsmooth
    ErrTolIntAccuracy = sim$para$eta/sim$para$rsmooth*1e-3 # => gives comparable accuracy to in-built simulator for Leapfrog integrator
    para = {}
    i = 0;
    i=i+1; para[i] = sprintf('InitCondFile                                      %s',filename.ics)
    i=i+1; para[i] = sprintf('OutputDir                                         %s',directory.output)
    i=i+1; para[i] = sprintf('SnapshotFileBase                                  snapshot')
    i=i+1; para[i] = sprintf('OutputListFilename                                empty.txt')
    i=i+1; para[i] = sprintf('ICFormat                                          1')
    i=i+1; para[i] = sprintf('SnapFormat                                        1')
    i=i+1; para[i] = sprintf('TimeLimitCPU                                      86400')
    i=i+1; para[i] = sprintf('CpuTimeBetRestartFile                             86400')
    i=i+1; para[i] = sprintf('MaxMemSize                                        2000')
    i=i+1; para[i] = sprintf('TimeBegin                                         0')
    i=i+1; para[i] = sprintf('TimeMax                                           %.12e',sim$para$t.max)
    i=i+1; para[i] = sprintf('ComovingIntegrationOn                             0')
    i=i+1; para[i] = sprintf('Omega0                                            0')
    i=i+1; para[i] = sprintf('OmegaLambda                                       0')
    i=i+1; para[i] = sprintf('OmegaBaryon                                       0')
    i=i+1; para[i] = sprintf('HubbleParam                                       1')
    i=i+1; para[i] = sprintf('Hubble                                            0')
    i=i+1; para[i] = sprintf('BoxSize                                           %.12e',sim$para$box.size)
    i=i+1; para[i] = sprintf('OutputListOn                                      0')
    i=i+1; para[i] = sprintf('TimeBetSnapshot                                   %.12e',sim$para$dt.out)
    i=i+1; para[i] = sprintf('TimeOfFirstSnapshot                               0')
    i=i+1; para[i] = sprintf('TimeBetStatistics                                 %.12e',sim$para$dt.out)
    i=i+1; para[i] = sprintf('NumFilesPerSnapshot                               1')
    i=i+1; para[i] = sprintf('MaxFilesWithConcurrentIO                          1')
    i=i+1; para[i] = sprintf('ErrTolIntAccuracy                                 %.12e',ErrTolIntAccuracy)
    i=i+1; para[i] = sprintf('CourantFac                                        0') # only used for SPH
    i=i+1; para[i] = sprintf('MaxSizeTimestep                                   %.12e',min(sim$para$dt.out,sim$para$dt.max))
    i=i+1; para[i] = sprintf('MinSizeTimestep                                   %.12e',0)
    i=i+1; para[i] = sprintf('TypeOfOpeningCriterion                            1')
    i=i+1; para[i] = sprintf('ErrTolTheta                                       0') # only used for tree scheme
    i=i+1; para[i] = sprintf('ErrTolThetaMax                                    0') # only used for tree scheme
    i=i+1; para[i] = sprintf('ErrTolForceAcc                                    0') # only used for tree scheme
    i=i+1; para[i] = sprintf('TopNodeFactor                                     2.5')
    i=i+1; para[i] = sprintf('ActivePartFracForNewDomainDecomp                  0.01')
    i=i+1; para[i] = sprintf('DesNumNgb                                         64')
    i=i+1; para[i] = sprintf('MaxNumNgbDeviation                                1')
    i=i+1; para[i] = sprintf('UnitLength_in_cm                                  1')
    i=i+1; para[i] = sprintf('UnitMass_in_g                                     1')
    i=i+1; para[i] = sprintf('UnitVelocity_in_cm_per_s                          1')
    i=i+1; para[i] = sprintf('GravityConstantInternal                           1')
    i=i+1; para[i] = sprintf('SofteningComovingClass0                           %.12e',smoothing)
    i=i+1; para[i] = sprintf('SofteningMaxPhysClass0                            %.12e',smoothing)
    i=i+1; para[i] = sprintf('SofteningClassOfPartType0                         0')
    i=i+1; para[i] = sprintf('SofteningComovingClass1                           %.12e',smoothing)
    i=i+1; para[i] = sprintf('SofteningMaxPhysClass1                            %.12e',smoothing)
    i=i+1; para[i] = sprintf('SofteningClassOfPartType1                         0')
    i=i+1; para[i] = sprintf('ArtBulkViscConst                                  1')
    i=i+1; para[i] = sprintf('MinEgySpec                                        0')
    i=i+1; para[i] = sprintf('InitGasTemp                                       0')
    write.table(para, file=filename.para, col.names = FALSE, row.names = FALSE, quote = FALSE)

    # run code
    if (is.null(sim$code$gadget.np)) sim$code$gadget.np = 1
    if (sim$code$gadget.np<1) stop('gadget.np must be a positive integer.')
    cmd = sprintf('mpirun -np %d %s %s', sim$code$gadget.np, sim$code$file, filename.para)
    error.code = system(cmd,intern=FALSE,ignore.stderr=TRUE,ignore.stdout=!verbose)
    if (error.code>0) stop('Could not execute the gadget code. Try recompiling gadget.')

    # load data
    n = length(sim$ics$m)
    n.snapshots = 0
    while(file.exists(sprintf('%s/snapshot_%03d',directory.output,n.snapshots))) n.snapshots=n.snapshots+1
    if (n.snapshots<2) stop('not enough GADGET-4 snapshots found')
    t.out = seq(0,by=sim$para$dt.out,length=n.snapshots)
    x.out = v.out = array(NA,c(n.snapshots,n,3))
    for (i in seq(n.snapshots)) {
      fn = sprintf('%s/snapshot_%03d',directory.output,i-1)
      dat = .gadget_read(fn)
      x.out[i,,] = dat$x
      v.out[i,,] = dat$v
    }

    # complete output
    sim$output = list(t = t.out, x = x.out, v = v.out,
                      n.snapshots = n.snapshots, n.iterations = NA,
                      n.acceleration.evaluations = NA)

  }


  # In-built simulator with C++ acceleration ************************************************************************

  else if (sim$code$name=='R') {

    if (sim$para$integrator=='euler') {

      iteration = function(dt) {
        .evaluate.accelerations()
        global$x = global$x+global$v*dt+0.5*global$a*dt^2
        global$v = global$v+global$a*dt
      }

    } else if (sim$para$integrator=='leapfrog') {

      iteration = function(dt) {
        global$v = global$v+global$a*dt/2
        global$x = global$x+global$v*dt
        .evaluate.accelerations()
        global$v = global$v+global$a*dt/2
      }

    } else if (sim$para$integrator%in%c('yoshida','yoshida4')) {

      .w1 = 1/(2-2^(1/3))
      .w0 = -2^(1/3)*.w1
      .c.yoshida = c(.w1/2,(.w0+.w1)/2,(.w0+.w1)/2,.w1/2)
      .d.yoshida = c(.w1,.w0,.w1)
      iteration = function(dt) {
        for (i in seq(3)) {
          global$x = global$x+.c.yoshida[i]*global$v*dt
          .evaluate.accelerations()
          global$v = global$v+.d.yoshida[i]*global$a*dt
        }
        global$x = global$x+.c.yoshida[4]*global$v*dt
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
          global$x = global$x+.c.yoshida[i]*global$v*dt
          .evaluate.accelerations()
          global$v = global$v+.d.yoshida[i]*global$a*dt
        }
        global$x = global$x+.c.yoshida[10]*global$v*dt
      }

    } else {

      stop('unknown integrator for this code')

    }

    # acceleration function
    .evaluate.accelerations = function() {
      global$n.acceleration.evaluations = global$n.acceleration.evaluations+1
      if (is.null(sim$para$afield)) {
        global$a[,] = 0
      } else {
        global$a = sim$para$afield(global$x,global$t)
      }
      f = accelerations(global$m,global$x,global$v,global$a,sim$para$G,rsmoothsqr,sim$para$box.size)
      global$a[,] = f$a
      global$dt.var = f$dtvar
    }

    # create a global environment (these are variables that can potentially be modified inside functions)
    global = new.env()
    global$m = sim$ics$m
    global$x = sim$ics$x
    global$v = sim$ics$v
    n = length(global$m)
    global$a = array(0,c(n,3))

    # initialize output variables
    n.out = ceiling(sim$para$t.max/sim$para$dt.out)+2
    global$t.out = rep(NA,n.out)
    if (sim$para$include.bg) {
      n.save = n
      save.list = 1:n
    } else {
      n.save = sum(sim$ics$m>=0)
      save.list = sim$ics$m>=0
    }
    if (n.save==0) stop('No particles in simulation output. Consider setting include.bg=TRUE.')
    global$x.out = global$v.out = array(NA,c(n.out,n.save,3))
    global$i.out = 0

    # initialize first iteration
    global$t = 0
    n.iterations = 0
    global$n.acceleration.evaluations = 0
    t.next = sim$para$dt.out # time of next output
    rsmoothsqr = sim$para$rsmooth^2
    .evaluate.accelerations()

    # save first snapshot
    .save.snapshot = function() {
      global$i.out = global$i.out+1
      global$t.out[global$i.out] = global$t
      global$x.out[global$i.out,,] = global$x[save.list,]
      global$v.out[global$i.out,,] = global$v[save.list,]
    }
    .save.snapshot()

    # iterate
    while (global$t < sim$para$t.max) {
      dt = min(sim$para$dt.max, sim$para$t.max-global$t, t.next-global$t, max(sim$para$dt.min, sim$para$eta*global$dt.var))
      iteration(dt)
      global$t = global$t+dt
      if (sim$para$box.size>0) global$x = global$x%%sim$para$box.size
      n.iterations = n.iterations+1
      if (global$t>=t.next | global$t>=sim$para$t.max) {
        .save.snapshot()
        t.next = global$i.out*sim$para$dt.out
      }
    }

    # finalise output
    sim$output = list(t = global$t.out[1:global$i.out], x = global$x.out[1:global$i.out,,], v = global$v.out[1:global$i.out,,],
                      n.snapshots = global$i.out, n.iterations = n.iterations, n.acceleration.evaluations = global$n.acceleration.evaluations)

  } else {

    stop(sprintf('unknown simulation code: %s',sim$code$name))

  }

  # turn list into a simulation class
  class(sim) = 'simulation'

  return(sim)

}

.gadget_write = function (ics, file, tinitial=0) {

  n = length(ics$m)

  data = file(file, "wb")

  # write header
  writeBin(as.integer(256), data, size = 4)
  writeBin(as.integer(c(0,n,0,0,0,0)), data, size = 4)
  writeBin(rep(0,6), data, size = 8)
  writeBin(as.numeric(tinitial), data, size = 8)
  writeBin(as.numeric(0), data, size = 8)
  writeBin(as.integer(rep(0,10)), data, size = 4)
  writeBin(as.numeric(rep(0,4)), data, size = 8)
  writeBin(as.integer(rep(0, 24)), data, size = 4)
  writeBin(as.integer(256), data, size = 4)

  # write positions
  writeBin(as.integer(n * 3 * 4), data, size = 4)
  writeBin(as.numeric(t(ics$x)), data, size = 4)
  writeBin(as.integer(n * 3 * 4), data, size = 4)

  # write velocities
  writeBin(as.integer(n * 3 * 4), data, size = 4)
  writeBin(as.numeric(t(ics$v)), data, size = 4)
  writeBin(as.integer(n * 3 * 4), data, size = 4)

  # write IDs
  writeBin(as.integer(n * 4), data, size = 4)
  writeBin(as.integer(seq(n)), data, size = 4)
  writeBin(as.integer(n * 4), data, size = 4)

  # write masses
  writeBin(as.integer(n * 4), data, size = 4)
  writeBin(as.numeric(ics$m), data, size = 4)
  writeBin(as.integer(n * 4), data, size = 4)

  close(data)
}

.gadget_read = function (file) {

  data = file(file, "rb")

  # read header
  block.start = readBin(data, "integer", n = 1, size=4)
  readBin(data, "integer", n = 1, size=4)
  n = readBin(data, "integer", n = 1, size=4)
  readBin(data, "integer", n = 62, size=4)
  block.end = readBin(data, "integer", n = 1, size=4)
  if (block.start!=block.end) {
    close(data)
    stop(paste0('incorrect block in file',file))
  }

  # read positions
  block.start = readBin(data, "integer", n = 1, size=4)
  x = t(array(readBin(data, "numeric", n = 3*n, size=4),c(3,n)))
  block.end = readBin(data, "integer", n = 1, size=4)
  if (block.start!=block.end) {
    close(data)
    stop(paste0('incorrect block in file',file))
  }

  # read velocities
  block.start = readBin(data, "integer", n = 1, size=4)
  v = t(array(readBin(data, "numeric", n = 3*n, size=4),c(3,n)))
  block.end = readBin(data, "integer", n = 1, size=4)
  if (block.start!=block.end) {
    close(data)
    stop(paste0('incorrect block in file',file))
  }

  # read IDs
  block.start = readBin(data, "integer", n = 1, size=4)
  id = readBin(data, "integer", n = n, size=4)
  block.end = readBin(data, "integer", n = 1, size=4)
  if (block.start!=block.end) {
    close(data)
    stop(paste0('incorrect block in file',file))
  }

  close(data)

  # sort particles by increasing id
  i = sort.int(id, index.return = T)$ix
  x = x[i,]
  v = v[i,]

  return(list(x=x, v=v))
}
