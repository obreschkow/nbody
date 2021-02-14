#' Show one or seveal 3D point sets in four projections
#'
#' @importFrom cooltools nplot runif3 rotation3
#' @importFrom grDevices pdf dev.off graphics.off
#' @importFrom graphics lines
#'
#' @description Produces a 2-by-2 panel figure, where each panel is made of an \code{\link{sphview}} image in a different rotation.
#'
#' @param x n-by-3 matrix, representing the 3D-coordinates of n particles.
#' @param rotations 4-vector or 4-element list specifying four different values/vectors of the argument \code{rotation} in \code{\link{sphview}}. These different rotations are shown in the bottom left, bottom right, top left, and top right panel, respectively, in this order.
#' @param screen logical flag specifying whether the images is displayed on the screen.
#' @param pdffile optional pdf-filename to save the image as pdf-file.
#' @param title Text to be added to the figure.
#' @param scale logical flag, specifying if a length scale is shown
#' @param ... additional parameters for \code{\link{sphview}}.
#'
#' @seealso \code{\link{sphview}}
#'
#' @author Danail Obreschkow
#'
#' @examples
#' # Example of a triaxial ellipsoid shown in four projections
#' x = t(t(cooltools::runif3(1e4))*c(3,2,1))%*%cooltools::rotation3(c(1,1,1))
#' sphview4(x)
#'
#' @export

sphview4 = function(x, rotations=c(2,3,4,1), screen = TRUE, pdffile = NULL, title = NULL, scale = TRUE, ...) {

  for (mode in seq(2)) {

    make = FALSE
    if (mode==1 & screen) {
      make = TRUE
    }
    if (mode==2 & !is.null(pdffile)) {
      make = TRUE
      grDevices::pdf(pdffile,width=7,height=7)
    }

    if (make) {

      if (mode==1) grDevices::graphics.off()
      cooltools::nplot(xlim=c(0,1), ylim=c(0,1), pty='s', mar=c(0,0,0,0))
      p = par()$plt

      ix = iy = 0

      for (i in seq(4)) {

        par(fig=c(p[1]+(p[2]-p[1])*ix/2,p[1]+(p[2]-p[1])*(ix/2+0.5),p[3]+(p[4]-p[3])*iy/2,p[3]+(p[4]-p[3])*(iy/2+0.5)),
            new=TRUE, mar=c(0,0,0,0) )

        if (i==1) {
          s = sphview(x, rotation = rotations[[i]], pngfile = NULL, pdffile = NULL, screen = TRUE,
                      scale=ifelse(i==2,scale,FALSE), ...)
        } else {
          s = sphview(x, rotation = rotations[[i]], pngfile = NULL, pdffile = NULL, screen = TRUE,
                      scale=ifelse(i==2,scale,FALSE), width=diff(s$xlim), ...)
        }
        ix = (ix+1)%%2
        if (ix==0) iy = iy+1

      }

      # lines between panels
      par(fig=p, new=TRUE, mar=c(0,0,0,0))
      cooltools::nplot(xlim=c(0,1), ylim=c(0,1))
      lines(c(1,1)/2,c(0,1),col='grey')
      lines(c(0,1),c(1,1)/2,col='grey')

      # title
      if (!is.null(title)) text(0.06, 1.94, title, pos=4, col='white', offset=-0.4)

      if (mode==2) grDevices::dev.off()

    }

  }

}
