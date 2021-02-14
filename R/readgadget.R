#' Read Gadget Data
#'
#' @importFrom rhdf5 h5read
#' @importFrom snapshot snapread
#'
#' @description Reads astrophysical N-body data output by the Gadget code
#'
#' @param file filename of file to load
#' @param type character specifying the data format. Must be either of: \code{bin} for binary format, \code{hdf} for HDF5 format or \code{auto} to automatically determine the format from file extension.
#'
#' @return Returns a list containing the particle data. The format of the list depends on the input format. Use \code{View(readgadget(...))} to see the output structure.
#'
#' @author Danail Obreschkow
#'
#' @export

readgadget = function(file, type='auto') {

  # check file
  if (!file.exists(file)) stop(paste0('File not found: ',file))
  if (file.access(file,4)!=0) stop(paste0('No permission to read: ',file))

  # handle type
  if (type=='binary') type = 'bin'
  if (type=='hdf5' | type=='h5') type = 'hdf'
  if (type=='auto') {

    # get file extension
    ext = strsplit(basename(file), split="\\.")[[1]]
    ext = ext[-1]

    # check file extension
    if (length(ext)>0) {
      if (ext%in%c('hdf','hdf5','h5')) {
        type = 'hdf'
      } else {
        type = 'bin'
      }
    } else {
      type = 'bin'
    }

  }

  # load file
  if (type=='bin') {
    return(snapshot::snapread(file))
  } else if (type=='hdf') {
    return(rhdf5::h5read(file,'/'))
  } else {
    stop('Unknown "type".')
  }

}
