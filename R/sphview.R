#' Visualise one or several 3D point sets
#'
#' @importFrom cooltools nplot rasterflip lim griddata2 kde2 quadrupole rotation3 vectornorm
#' @importFrom png writePNG
#' @importFrom EBImage gblur
#' @importFrom grDevices pdf dev.off col2rgb
#' @importFrom graphics axis lines par rasterImage rect text
#'
#' @description Produces a raster image visualising particle positions in 3D N-body/SPH data.
#'
#' @param x n-by-3 matrix, representing the 3D-coordinates of n particles.
#' @param species optional argument. If \code{x} is a n-by-3 matrix, \code{species} is a n-element vector giving the integer species indices, which will specify the color of the particles; if not given, all particles will be assumed to be of species=1. If \code{x} is a \code{snapshot} list, \code{species} is a single integer or a vector of integer specifying the species to select; if not given, all species are plotted.
#' @param value optional n-element vector specifying a particle property (e.g. temperatures) to be shown in color. This overwrites the colors set by \code{species}.
#' @param valrange 2-vector specifying the range of values covered by the colors in the vector \code{colorscale}.
#' @param fix.brightness logical flag specifying whether the brightness of value-plots should scale with density.
#' @param weight optional n-element vector, specifying the weight of different particles.
#' @param colspecies vector of colors, specifying the color of the species given in the optional \code{species} vector.
#' @param colscale vector of colors used to shown the particle property given by the optional vectors \code{value}.
#' @param center optional 3-vector specifying the coordinate at the center of the plot. The default is the geometric center (= center of mass, if all particle masses are equal), if \code{x} is a matrix; and the box center, if \code{x} is a Gadget snapshot list.
#' @param width optional horizontal range of the image in the length units of \code{x}. The default corresponds to the full range of particle positions, if \code{x} is a matrix; and to the box size, if \code{x} is a Gadget snapshot list.
#' @param aspect aspect ratio (= width/height).
#' @param thickness optional thickness in length units of \code{x} of the slice to be plotted. If not given, all particles are shown.
#' @param subsampling fraction of particles to be used (if 1, all particles are used, if <1 a random subsample is drawn)
#' @param screen logical flag specifying whether the images is displayed on the screen.
#' @param pngfile optional png-filename to save the image as raster image.
#' @param pdffile optional pdf-filename to save the image as pdf-file.
#' @param rotation either an integer (1-6), a 3-vector or a 3-by-3 matrix, specifying a rotation of the 3D particle positions. In case of an integer: 1=(x,y)-plane, 2=(y,z)-plane, 3=(x,y)-plane, 4=(qmax,qmin)-plane, 5=(qmax,qmid)-plane, 6=(qmid,qmin)-plane, where qmax/qmid/qmin are the eigenvectors of the particle-quadrupole, associated with the maximum/middle/minimum eigenvalues, respectively. If \code{rotation} is a vector, its direction specifies the rotation axis and its norm the rotation angle in radians in the positive geometric sense.
#' @param ngrid optional number of grid cells per side in the output image. If the image is not square ngrid is interpreted as the geometric mean between the horizontal and the vertical number of pixels, such that the total number of pixels remains about ngrid^2.
#' @param title Text to be added to the figure.
#' @param title.origin optional 2-vector specifying the position of the title
#' @param lum overall luminosity scaling factor (default 1).
#' @param shadows differential luminosity scaling factor for darker regions (default 1).
#' @param freeze.lum logical flag (default FALSE). If set to TRUE, the luminosity will not be adjusted as a function of the visible particles. This is normally needed to avoid flickering in animations.
#' @param sigma smoothing scale in simulation units. The default is NULL, meaning that a value is picked automatically. If a vector is provided, the different elements will be interpreted as corresponding to different species.
#' @param kde logical flag, specifying whether the particles are smoothed using an adaptive kernel density estimator. If a vector is provided, the different elements will be interpreted as corresponding to different species.
#' @param arrows logical flag, specifying if axis-arrows are drawn
#' @param arrow.origin optional 2-vector specifying the bottom left corner of the arrows
#' @param arrow.length optional number specifying the length of the arrows
#' @param arrow.lwd line width of arrows
#' @param scale logical flag, specifying if a length scale is shown
#' @param scale.origin optional 2-vector specifying the right end of the length scale
#' @param scale.length optional number specifying the length of the length scale (is always rounded to one significant digit)
#' @param scale.lwd line width of length scale
#' @param length.unit character string of length unit (e.g. "m" or "kpc")
#' @param xlab label on horizontal arrow (only shown if \code{arrows=TRUE}).
#' @param ylab label on vertical arrow (only shown if \code{arrows=TRUE}).
#' @param cex text size
#' @param text.offset optional extra offset for the arrow and scale labels
#'
#' @return Returns a list of items
#' \item{rgb}{ngrid-by-ngrid-by-3 array of containig the rgb colour data between 0 and 1}
#' \item{radius}{positive real number specifying the radius of the image along the horizontal and vertical axes, in the units of the particle positions}
#' \item{rotationmatrix}{rotation matrix applied to the input data, i.e. the projection was x = (x\%*\%rot)[,1:2]}
#'
#' @author Danail Obreschkow
#'
#' @examples
#' # Example of 2x1e4 particles in a triaxial ellipsoid
#' x.red = t(t(cooltools::runif3(1e4))*c(3,2,1))
#' x.blue = abs(t(t(cooltools::runif3(1e4,r=c(0.5,1)))*c(3,2,1)))
#' sphview(rbind(x.red,x.blue), species=c(rep(1,1e4),rep(2,1e4)), center=0)
#'
#' @export

sphview = function(x, species = NULL, value = NULL, valrange = NULL, fix.brightness = TRUE, weight = NULL,
                   colspecies = c('#ff0010', '#0515ff', 'green', 'orange', 'yellow', 'purple'),
                   colscale = cooltools::spectrumcolors(1000),
                   center = NULL, width = NULL, aspect = 1, thickness = NULL, subsampling = 1, screen = TRUE, pngfile = NULL, pdffile = NULL,
                   rotation = 1, ngrid = 300, title = NULL, title.origin = NULL,
                   lum = 1, shadows = 1, freeze.lum = FALSE, sigma = NULL, kde = TRUE,
                   arrows = TRUE, arrow.origin = NULL, arrow.length = NULL, arrow.lwd = 1.5,
                   scale = TRUE, scale.origin = NULL, scale.length = NULL, scale.lwd = 1.5,
                   length.unit = '', xlab = NULL, ylab = NULL, cex=1, text.offset = 0) {


  # handle x
  if (is.array(x)) {
    if (length(dim(x))==2) {
      n = dim(x)[1]
      if (dim(x)[2]==2) x = cbind(x,rep(0,n))
      if (dim(x)[2]!=3) stop('x must be an n-by-3 dimensional array')
    } else {
      stop('x must be an n-by-3 dimensional array')
    }
  } else {
    stop('x must be an n-by-3 dimensional array')
  }

  # handle species
  if (!is.null(species)) {
    if (length(species)!=n) stop('species must be a vector of same length as number of rows in x.')
    if (min(species)<1) stop('species must contain only positive integers.')
    if (max(species)>length(colspecies)) stop('not enough colours given to match the maximum species value.')
  }

  # handle values
  if (!is.null(value)) {
    if (length(value)!=n) stop('value must be a vector of same length as number of rows in x.')
    if (is.null(valrange)) valrange = range(value)+c(-1e-5,1e-5)*(mean(abs(value))+1)
  }

  # handle weights
  if (!is.null(weight)) {
    if (length(weight)!=n) stop('weight must be a vector of same length as number of rows in x.')
  }

  # subsampling
  if (subsampling<=0 | subsampling>1) stop('subsampling must be a value in the interval (0,1].')
  if (subsampling<1) {
    n.new = max(1,round(n*subsampling))
    id = sample(n,n.new)
    x = x[id,]
    if (!is.null(species)) species = species[id]
    if (!is.null(value)) value = value[id]
    if (!is.null(weight)) weight = weight[id]
    n = n.new
  }

  # center points
  if (is.null(center)) {
    center = apply(x,2,mean) # geometric centre
  } else {
    if (length(center)==2) center=c(center,0)
  }
  x = sweep(x, 2, center)

  # rotation matrix
  if (length(rotation)==3) {
    rot = t(cooltools::rotation3(rotation))
    if (is.null(xlab)) xlab = expression(e[1])
    if (is.null(ylab)) ylab = expression(e[2])
  } else if (length(rotation)==1) {
    if (rotation>3) e = eigen(cooltools::quadrupole(x))$vectors
    if (rotation==1) {
      rot = diag(3)
      if (is.null(xlab)) xlab = 'x'
      if (is.null(ylab)) ylab = 'y'
    } else if (rotation==2) {
      rot = rbind(c(0,1,0),c(0,0,1),c(1,0,0))
      if (is.null(xlab)) xlab = 'y'
      if (is.null(ylab)) ylab = 'z'
    } else if (rotation==3) {
      rot = rbind(c(0,0,1),c(1,0,0),c(0,1,0))
      if (is.null(xlab)) xlab = 'z'
      if (is.null(ylab)) ylab = 'x'
    } else if (rotation==4) {
      rot = t(e[,c(1,3,2)])
      if (is.null(xlab)) xlab = expression(lambda ['max'])
      if (is.null(ylab)) ylab = expression('  '*lambda ['min'])
    } else if (rotation==5) {
      rot = t(e[,c(1,2,3)])
      if (is.null(xlab)) xlab = expression(lambda ['max'])
      if (is.null(ylab)) ylab = expression('  '*lambda ['mid'])
    } else if (rotation==6) {
      rot = t(e[,c(2,3,1)])
      if (is.null(xlab)) xlab = expression(lambda ['mid'])
      if (is.null(ylab)) ylab = expression('  '*lambda ['min'])
    } else {
      stop('rotation must be an integer 1,...,6, a real 3-vector, or a 3-by-3 matrix.')
    }
  } else if (length(rotation)==9) {
    rot = matrix(rotation,3,3)
  } else {
    stop('rotation must be an integer 1,...,6, a real 3-vector, or a 3-by-3 matrix.')
  }

  # rotate coordinates (active rotation)
  x = t(rot%*%t(x))

  # select slice
  if (!is.null(thickness)) {
    sel = abs(x[,3])<=thickness/2
    n = sum(sel)
    if (n==0) stop('No particles to plot. Consider increasing "thickness".')
    x = x[sel,]
    if (!is.null(species)) species = species[sel]
    if (!is.null(value)) value = value[sel]
    if (!is.null(weight)) weight = weight[sel]
  }

  # determine plotting limits
  if (is.null(width)) width = 2*sqrt(2)*max(vectornorm(x))/sqrt(1+1/aspect^2)
  height = width/aspect
  xlim = c(-1,1)*width/2
  ylim = c(-1,1)*height/2
  mean.length = sqrt(width*height)

  # project coordinates
  x = x[,1:2]
  n.plot = sum(x[,1]>=xlim[1] & x[,1]<=xlim[2] & x[,2]>=ylim[1] & x[,2]<=ylim[2])

  # determine arrows, scale, titles
  if (is.null(arrow.origin)) arrow.origin = c(xlim[1],ylim[1])+0.05*mean.length
  if (is.null(arrow.length)) arrow.length = 0.1*mean.length
  if (is.null(scale.origin)) scale.origin = c(xlim[2]-0.05*mean.length,ylim[1]+0.05*mean.length)
  if (is.null(scale.length)) scale.length = 0.1*mean.length
  scale.length = signif(scale.length,1)
  if (is.null(title.origin)) title.origin = c(xlim[1]+0.06*mean.length,ylim[2]-0.06*mean.length)

  # prepare grid
  nx = round(ngrid*width/mean.length)
  ny = round(ngrid*height/mean.length)
  if (max(nx,ny)>8000) stop('Not more than 8000 pixels allowed on each side.')
  dx = width/nx # pixel size in simulation length units
  if (is.null(sigma)) {
    smoothing = 0.008*ngrid
  } else {
    smoothing = sigma/dx # smoothing length in pixel
  }
  img = array(0,c(nx,ny,3))

  # prepare rastering
  if (is.null(value) & is.null(species)) {
    type = 0
    s = 1
    n.loop = 1
    col = col2rgb(colspecies[1])/255
  } else if (is.null(value) & !is.null(species)) {
    type = 1
    s = unique(species)
    n.loop = length(s)
    col = t(col2rgb(colspecies)/255)
  } else {
    type = 2
    n.loop = 2
    col = t(col2rgb(colscale)/255)
  }

  if (length(smoothing)==1) smoothing=rep(smoothing,n.loop)
  if (length(kde)==1) kde=rep(kde,n.loop)

  n.eff = 0
  w = weight
  sel = seq(n)

  for (i in seq(n.loop)) {

    if (type == 1) {
      sel = species==s[i]
    } else if (type == 2) {
      if (i==1) {
        w = weight
      } else {
        w = value
        if (!is.null(weight)) w = w*weight
      }
    }

    if (smoothing[i]==0) {
      g = cooltools::griddata2(x[sel,1],x[sel,2],w=w[sel],xlim=xlim,ylim=ylim,n=c(nx,ny))
      if (is.null(g$m)) g$m = g$n
    } else {
      if (kde[i]) {
        g = cooltools::kde2(x[sel,1],x[sel,2],w=w[sel],xlim=xlim,ylim=ylim,n=c(nx,ny),s = smoothing[i]/8, sd.max=smoothing[i]*2)
        g$m = g$d
      } else {
        g = cooltools::griddata2(x[sel,1],x[sel,2],w=w[sel],xlim=xlim,ylim=ylim,n=c(nx,ny))
        if (is.null(g$m)) g$m = g$n
        g$m = EBImage::gblur(g$m,smoothing[i])
      }
    }

    if (type==0) {
      for (k in seq(3)) {
        img[,,k] = img[,,k]+col[k]*g$m
      }
      n.eff = n.eff+sum(g$m)
    } else if (type==1) {
      for (k in seq(3)) {
        img[,,k] = img[,,k]+col[s[i],k]*g$m
      }
      n.eff = n.eff+sum(g$m)
    } else if (type==2) {
      if (i==1) {
        img[,,1] = g$m
        n.eff = n.eff+sum(g$m)
      } else {
        meanval = g$m/img[,,1]
        meanval[!is.finite(meanval)] = NA
      }
    }

  }

  # rescale brightness
  img = atan(img*ngrid^2/ifelse(freeze.lum,n,n.eff)*(0.4*lum)*ifelse(type==2,0.3,0.5))/pi*2
  f = 10^max(0,2.1*shadows)
  img = log10(f*img+1)/log10(f+1)
  img = cooltools::lim(img) # necessary because gblur sometimes produces very slightly negative values
  if (type==2) {
    if (fix.brightness) {
      l = array(0.8,c(nx,ny))
    } else {
      l = img[,,1]
    }
    colid = pmin(length(colscale),pmax(1,round((meanval-valrange[1])/diff(valrange)*length(colscale)+0.5)))
    colid[is.na(colid)] = length(colscale)+1
    col = rbind(col,c(0,0,0))
    for (k in seq(3)) {
      img[,,k] = l*array(col[colid,k],c(nx,ny))
    }
  }

  # save raster image as png
  if (!is.null(pngfile)) {
    png::writePNG(cooltools::rasterflip(img),pngfile)
  }

  # show on screen and save as pdf
  for (mode in seq(2)) {

    make = FALSE

    if (mode==1 & screen) {
      make = TRUE
    }

    if (mode==2 & !is.null(pdffile)) {
      make = TRUE
      grDevices::pdf(pdffile,width=7*width/mean.length,height=7*height/mean.length)
    }

    if (make) {

      # initialize plot
      cooltools::nplot(xlim=xlim, ylim=ylim, asp=1)

      # plot raster
      rasterImage(cooltools::rasterflip(img),xlim[1],ylim[1],xlim[2],ylim[2])

      # arrows
      if (arrows) {
        arrows(arrow.origin[1],arrow.origin[2],arrow.origin[1]+arrow.length,arrow.origin[2],col='white',length = 0.1,angle=20,lwd=arrow.lwd)
        arrows(arrow.origin[1],arrow.origin[2],arrow.origin[1],arrow.origin[2]+arrow.length,col='white',length = 0.1,angle=20,lwd=arrow.lwd)
        if (!is.null(xlab)) {
          text(arrow.origin[1]+arrow.length+text.offset,arrow.origin[2],xlab,pos=4,col='white',cex=cex)
        }
        if (!is.null(ylab)) {
          text(arrow.origin[1],arrow.origin[2]+arrow.length+text.offset,ylab,pos=3,col='white',cex=cex)
        }
      }

      # length scale
      if (scale) {
        text(scale.origin[1]-scale.length/2,scale.origin[2]+text.offset,sprintf('%s %s',signif(scale.length,1),length.unit),pos=3,col='white',cex=cex)
        lines(scale.origin[1]-c(0,scale.length),rep(scale.origin[2],2),col='white',lwd=scale.lwd)
      }

      # title
      if (!is.null(title)) text(title.origin[1], title.origin[2], title, pos=4, col='white', offset=-0.4)

      if (mode==2) {dev.off()}

    }

  }

  # return results
  invisible(list(rgb = img, xlim = xlim, ylim = ylim, center=center, rotationmatrix = rot))

}
