#' @title Mechanical energy of an N-body system
#'
#' @description Computes the instantaneous potential and kinetic energies of all particles in an N-body system. Here, the potential energy of a particle i means the potential energy it has with all other particles (sum_j -G*m[i]*[j]/rij). Hence the total potential energy of the system is half the sum of the individual potential energies.
#'
#' @param m N-vector with the masses of the N particles. Negative masses are treated as positive masses of same magnitude, since negative masses normally represent positive background masses in the nbody package.
#' @param x N-by-3 matrix specifying the initial position in Cartesian coordinates
#' @param v N-by-3 matrix specifying the initial velocities
#' @param G gravitational constant. The default is the measured value in SI units.
#' @param rsmooth top-hat smoothing radius. If 0 given, no smoothing is assumed.
#' @param box.size scalar>=0. If 0, open boundary conditions are adopted. If >0, potential energies are computed assuming a cubic box of side length box.size with periodic boundary conditions. In this case, the cubic box is contained in the interval [0,box.size) in all three Cartesian coordinates. The separation between any two particles is always calculated along their shortest separation, which may cross 0-3 boundaries.
#' @param cpp logical flag. If TRUE (default), the computation is performed efficiently in C++.
#'
#' @return Returns a list with vector items \code{Ekin}, \code{Epot}, \code{Emec=Ekin+Epot}; and the associated total quantities \code{Ekin.tot}, \code{Epot.tot}, \code{Emec=Ekin+Epot.tot}.
#'
#' @author Danail Obreschkow
#'
#' @export

energy = function(m,x,v,rsmooth=0,box.size=0,G=6.67408e-11,cpp=TRUE) {

  m = abs(m) # to make sure that fixed masses

  if (dim(x)[2]!=3) stop('x must be a matrix with 3 columns.')
  if (dim(v)[2]!=3) stop('v must be a matrix with 3 columns.')
  if (dim(x)[1]!=dim(v)[1]) stop('x and v must have the same number of rows.')
  if (dim(x)[1]!=length(m)) stop('The length of m must be equal to the number of rows of x.')
  if (dim(v)[1]!=length(m)) stop('The length of m must be equal to the number of rows of v.')

  n = length(m)

  Ekin = 0.5*m*rowSums(v^2)
  Ekin.tot=sum(Ekin)

  if (cpp) {

    Epot = energypot(m,x,v,rsmooth,box.size)$Epot

  } else {

    Epot = rep(0,n)

    if (box.size==0) {

      for (i in seq(1,n-1)) {
        for (j in seq(i+1,n)) {
          r = sqrt(sum((x[i,]-x[j,])^2))
          if (r>rsmooth) {
            dEpot = -m[i]*m[j]/r
          } else {
            dEpot = -m[i]*m[j]*(3-r^2/rsmooth^2)/(2*rsmooth)
          }
          Epot[i] = Epot[i]+dEpot
          Epot[j] = Epot[j]+dEpot
        }
      }

    } else if (box.size>0) {

      h = box.size/2
      for (i in seq(1,n-1)) {
        for (j in seq(i+1,n)) {
          r = sqrt(sum(((x[i,]-x[j,]+h)%%box.size-h)^2))
          if (r>rsmooth) {
            dEpot = -m[i]*m[j]/r
          } else {
            dEpot = -m[i]*m[j]*(3-r^2/rsmooth^2)/(2*rsmooth)
          }
          Epot[i] = Epot[i]+dEpot
          Epot[j] = Epot[j]+dEpot
        }
      }

    } else {

      stop('Invalid box.size.')

    }

  }

  Epot = Epot*G
  Epot.tot = sum(Epot)/2

  return(list(Ekin=Ekin,Epot=Epot,Emec=Ekin+Epot,
              Ekin.tot=Ekin.tot,Epot.tot=Epot.tot,Emec.tot=Ekin.tot+Epot.tot))
}
