#' @title Visualize an N-body simulation
#'
#' @description Basic routine to visualise the result of an N-body simulation, projected onto a plane.
#'
#' @importFrom magicaxis magplot
#' @importFrom graphics lines par points
#' @importFrom grDevices col2rgb
#'
#' @param x is a simulation-object as produced by \code{run.simulation}
#' @param y deprecated argument included for consistency with generic \code{\link[graphics]{plot}} function
#' @param units length unit in SI units
#' @param index1 index of the dimension plotted on the x-axis
#' @param index2 index of the dimension plotted on the y-axis
#' @param xlim 2-vector specifying the plotting range along the x-axis
#' @param ylim 2-vector specifying the plotting range along the y-axis
#' @param center 3-vector specifying the plotting center in the specified units
#' @param cex point size
#' @param pch point type
#' @param title title of plot
#' @param asp aspect ratio of x and y axes
#' @param pty character specifying the type of plot region to be used; "s" generates a square plotting region and "m" generates the maximal plotting region.
#' @param col either (1) a single color, (2) a n-element vector of colors for each particle or (3) a function(n,...) producing n colors, e.g. 'rainbow'
#' @param alpha.orbits opacity (0...1) of orbital lines.
#' @param alpha.snapshots opacity (0...1) of snapshot points.
#' @param lwd line width of orbital lines.
#' @param show.orbits logical flag. If TRUE (default), the orbits are shown as straight lines between snapshots.
#' @param show.snapshots logical flag. If TRUE (default), points are shown for each snapshot.
#' @param show.ics logical flag. If TRUE (default), the initial positions are highlighted.
#' @param show.fcs logical flag. If TRUE (default), the final positions are highlighted.
#' @param ... additional parameters for \code{\link[graphics]{plot}}
#'
#' @author Danail Obreschkow
#'
#' @method plot simulation
#' @export
plot.simulation = function(x, y, units=1, index1=1, index2=2, xlim=NULL, ylim=NULL,
                           center=c(0,0,0), cex=0.3, pch=20, title='', asp=1, pty='m', col='black',
                           alpha.orbits=1, alpha.snapshots=1, lwd=1,
                           show.orbits=TRUE, show.snapshots=TRUE, show.ics=TRUE, show.fcs=TRUE, ...) {

  # input check
  if (is.null(x$output)) stop('It seems that this simulation has not yet been run, as its output list is missing.')
  if (length(dim(x$output$x))==3) {
    n = dim(x$output$x)[2] # number of particles
  } else if (length(dim(x$output$x))==2) {
    n = 1
    x$output$x = array(x$output$x,c(dim(x$output$x)[1],1,dim(x$output$x)[2]))
    x$output$v = array(x$output$v,c(dim(x$output$v)[1],1,dim(x$output$v)[2]))
  } else {
    stop('unknown data structure')
  }

  # handle colors
  if (is.function(col)) {
    col = col(n)
  } else if (length(col)==1) {
    col = rep(col,n)
  } else {
    if (length(col)!=n) stop('col must be a single color or a n-element vector of colors.')
  }

  # produce coordinates
  x = x$output$x/units
  for (i in seq(3)) x[,,i] = x[,,i]-center[i]

  # open empty plot
  if (is.null(xlim)) xlim = range(x[,,index1])
  if (is.null(ylim)) ylim = range(x[,,index2])
  par(pty=pty)
  magplot(NA,xlim=xlim,ylim=ylim,asp=asp,...)

  # draw trajectories
  m = dim(x)[1]
  for (i in seq(n)) {
    if (show.orbits) lines(x[,i,index1],x[,i,index2],col=.transparent(col[i],alpha.orbits),lwd=lwd)
    if (show.snapshots) points(x[,i,index1],x[,i,index2],pch=pch,cex=cex,col=.transparent(col[i],alpha.snapshots))
    if (show.ics) points(x[1,i,index1],x[1,i,index2],pch=1,cex=cex*3,col=col[i])
    if (show.fcs) points(x[m,i,index1],x[m,i,index2],pch=4,cex=cex*2.2,col=col[i])
  }

  # finalize plot
  title(title,cex.main=0.8,line=0)

}

.transparent = function(col,alpha=0.5) {
  for (i in seq(length(col))) {
    rgb = col2rgb(col[i],alpha=TRUE)/255
    col[i] = rgb(rgb[1],rgb[2],rgb[3],rgb[4]*alpha)
  }
  return(col)
}
