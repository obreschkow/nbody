#' Set a default external simulation code
#'
#' @param code structured list specifying the default external simulation code used when calling \code{\link{run.simulation}}. This list has exactly the same format as the sub-list `code` described in the documentation of \code{\link{run.simulation}}.
#'
#' @return Returns the current list `code`. If no such last has been set and `default.code()` is called without argument, an error is produced.
#'
#' @author Danail Obreschkow
#'
#' @export

default.code = function(code=NULL) {

  if (is.null(code)) {
    if (is.null(.nbody.env$code)) {
      stop('no default simulation code has been set, yet. call default.code with an argument.')
    } else {
      code = .nbody.env$code
    }
  } else {
    .nbody.env$code = code
  }

  invisible(.nbody.env$code)
}
