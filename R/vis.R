# generic routine to plot the orbits computed with `run.simulation`
plot.simulation = function(sim, units, cex=0.3, pch=20, title='', ...) {
  par(pty='s')
  magplot(NA,NA, ...)
  x.out = sim$output$x.out/units
  n.out = dim(x.out)[1]
  for (i in seq(dim(x.out)[3])) {
    lines(x.out[,1,i],x.out[,2,i])
    points(x.out[,1,i],x.out[,2,i],pch=pch,cex=cex)
    points(x.out[1,1,i],x.out[1,2,i],pch=1,cex=1,col='blue')
    points(x.out[n.out,1,i],x.out[n.out,2,i],pch=20,cex=0.5,col='blue')
  }
  title(title,cex.main=0.8,line=0)
}