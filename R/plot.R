#' @title Visualize an N-body simulation
#'
#' @description Basic routine to visualise the result of an N-body simulation, projected onto a plane.
#'
#' @importFrom magicaxis magplot
#' @importFrom graphics lines par points
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
#' @param ... additional parameters for \code{\link[graphics]{plot}}
#'
#' @author Danail Obreschkow
#'
#' @method plot simulation
#' @export
plot.simulation = function(x, y, units=1, index1=1, index2=2, xlim=NULL, ylim=NULL, center=c(0,0,0), cex=0.3, pch=20, title='', asp=1, pty='m', ...) {

  # input check
  if (is.null(x$output)) stop('It seems that this simulation has not yet been run, as its output list is missing.')

  # produce coordinates
  x = x$output$x/units
  for (i in seq(3)) x[,,i] = x[,,i]-center[i]

  # open empty plot
  par(pty=pty)
  if (is.null(xlim)) xlim = range(x[,,index1])
  if (is.null(ylim)) ylim = range(x[,,index2])
  magplot(NA,xlim=xlim,ylim=ylim,asp=asp,...)

  # draw trajectories
  n = dim(x)[1]
  for (i in seq(dim(x)[2])) {
    lines(x[,i,index1],x[,i,index2])
    points(x[,i,index1],x[,i,index2],pch=pch,cex=cex)
    points(x[1,i,index1],x[1,i,index2],pch=1,cex=1,col='blue')
    points(x[n,i,index1],x[n,i,index2],pch=20,cex=0.5,col='blue')
  }

  # finalize plot
  title(title,cex.main=0.8,line=0)

}
