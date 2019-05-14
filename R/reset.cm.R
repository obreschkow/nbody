#' @title Move center of mass to the origin
#'
#' @description Routine, designed to reset the center of mass (CM) of the initial conditions (ICs) of an N-body simulation. The CM position and velocity are both shifted to (0,0,0).
#'
#' @param sim list of m, x, v or list with a sublist "ics", made of m, x, v, where
#' \code{m} = N-vector with the masses of the N particles\cr
#' \code{x} = N-by-3 matrix specifying the initial position in cartesian coordinates\cr
#' \code{v} = N-by-3 matrix specifying the initial velocities\cr\cr
#'
#' @return Returns a structure of the same format as the input argument, but with re-centered positions and velocities.
#'
#' @author Danail Obreschkow
#'
#' @export

reset.cm = function(sim) {

  center = function(m,x,v) {

    if (dim(x)[2]!=3) stop('x must be a matrix with 3 columns.')
    if (dim(v)[2]!=3) stop('v must be a matrix with 3 columns.')
    if (dim(x)[1]!=dim(v)[1]) stop('x and v must have the same number of rows.')
    if (dim(x)[1]!=length(m)) stop('The length of m must be equal to the number of rows of x.')
    if (dim(v)[1]!=length(m)) stop('The length of m must be equal to the number of rows of v.')

    M = sum(m)
    if (M<=0) stop('Total mass must be positive.')

    x0 = colSums(x*m)/M # center of mass position
    v0 = colSums(v*m)/M # center of mass velocity
    x = t(t(x)-x0)
    v = t(t(v)-v0)

    return(list(m=m, x=x, v=v))

  }

  if (is.null(sim$ics)) {

    if (is.null(sim$m)) stop('"sim$m" must be given is "sim$ics" is not given.')
    if (is.null(sim$x)) stop('"sim$x" must be given is "sim$ics" is not given.')
    if (is.null(sim$v)) stop('"sim$v" must be given is "sim$ics" is not given.')
    return(center(sim$m,sim$x,sim$v))

  } else {

    if (is.null(sim$ics$m)) stop('"sim$ics$m" must be given if "sim$ics" is given.')
    if (is.null(sim$ics$x)) stop('"sim$ics$x" must be given if "sim$ics" is given.')
    if (is.null(sim$ics$v)) stop('"sim$ics$v" must be given if "sim$ics" is given.')
    sim$ics = center(sim$ics$m,sim$ics$x,sim$ics$v)
    return(sim)

  }

}
