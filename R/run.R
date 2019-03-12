# generic routine to run an N-body simulation
run.sim = function(sim) {
  
  # make global variables from ICs
  m <<- sim$ics$m
  x <<- sim$ics$x
  v <<- sim$ics$v
  n <<- length(m)
  
  # set iteration function
  if (is.null(sim$para$algorithm)) sim$para$algorithm = 'leapfrog'
  custom.iteration = iteration[[sim$para$algorithm]]
  
  # initialize output variables
  n.out = floor(sim$para$t.max/sim$para$dt.out)+2
  x.out = v.out = array(NA,c(n.out,3,n))
  i.out = 1
  x.out[1,,] = x
  v.out[1,,] = v
  t.out = sim$para$dt.out # time of next output
  
  # advanced settings (for future use)
  rsmoothsqr <<- sim$para$rsmooth^2 # only used for smoothed particles
  dt.var <<- NULL # only used for variable time-stepping
  
  # prepare first iteration
  t = 0
  n.iterations = 0
  evaluate.accelerations()
  
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

# Wrapper to catch errors and measure simulation time
run.simulation = function(sim) {
  return(tryCatch({
    t.start = Sys.time()
    sim = run.sim(sim)
    t.end = Sys.time()
    time.taken = as.double(t.end)-as.double(t.start)
    message(sprintf('Simulation successfully completed in %.2fs.',time.taken))
  }, warning = function(w) {
    message('Warning(s) produced in executing the simulation.')
    sim$output = NA
  }, error = function(e) {
    message('Error(s) produced in executing the simulation.')
    sim$output = NA
  }, finally = {
    return(sim)
  }))
}

# generic acceleration evaluation in C
code = '
for (int i = 0; i<(*n*3); i++) a[i] = 0;
double dminsqr[*n];
for (int i = 0; i<(*n); i++) dminsqr[i] = 1e99;
for (int i = 0; i<(*n-1); i++) {
int i3 = 3*i;
for (int j = i+1; j<(*n); j++) {
if ((m[i]>0)||(m[j]>0)) {
int j3 = 3*j;
double dx = x[i3]-x[j3];
double dy = x[i3+1]-x[j3+1];
double dz = x[i3+2]-x[j3+2];
double dsqr = fmax(*rsmoothsqr,dx*dx+dy*dy+dz*dz);
double q = *G/pow(dsqr,1.5);
double qx = q*dx;
double qy = q*dy;
double qz = q*dz;
a[i3]-=m[j]*qx; a[i3+1]-=m[j]*qy; a[i3+2]-=m[j]*qz;
a[j3]+=m[i]*qx; a[j3+1]+=m[i]*qy; a[j3+2]+=m[i]*qz;
dminsqr[i] = fmin(dminsqr[i],dsqr);
dminsqr[j] = fmin(dminsqr[j],dsqr);
}
}
}
double zmin = 1e99;
for (int i = 0; i<(*n); i++) {
int i3 = 3*i;
double z = dminsqr[i]/(a[i3]*a[i3]+a[i3+1]*a[i3+1]+a[i3+2]*a[i3+2]);
if (z<zmin) zmin=z;
}
*dtvar = pow(zmin,0.25);
'

evaluate.accelerations.c = cfunction(signature(n='i',m='n',x='n',a='n',
                                               dtvar='n',rsmoothsqr='n',G='n'),code,convention=".C")

evaluate.accelerations = function() {
  a <<- array(0,c(3,n))
  f = evaluate.accelerations.c(n,m,x,a,0,rsmoothsqr,cste$G)
  a[,] <<- f$a
  dt.var <<- f$dtvar
}

# basic iteration
iteration = {}
iteration$euler = function(dt) {
  evaluate.accelerations()
  x <<- x+v*dt+0.5*a*dt^2
  v <<- v+a*dt
}

iteration$leapfrog = function(dt) {
  v <<- v+a*dt/2
  x <<- x+v*dt
  evaluate.accelerations()
  v <<- v+a*dt/2
}

w1 = 1/(2-2^(1/3))
w0 = -2^(1/3)*w1
c.yoshida = c(w1/2,(w0+w1)/2,(w0+w1)/2,w1/2)
d.yoshida = c(w1,w0,w1)
iteration$yoshida = function(dt) {
  for (i in seq(3)) {
    x <<- x+c.yoshida[i]*v*dt
    evaluate.accelerations()
    v <<- v+d.yoshida[i]*a*dt
  }
  x <<- x+c.yoshida[4]*v*dt
}