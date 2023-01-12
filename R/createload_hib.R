#' Create a database object
#' 
#' @description 
#' This function creates a [`R6`] object containing the loaded \code{PREDICTS} database. It
#' assumes that both a species occurrence and a site database is present. All files are loaded and
#' compiled in a single database that can be stored for future use.
#' 
#' This object is then used in other functions of this package to prepare or filter the
#' database and enable ecological impact analyses. 
#' @note 
#' Other or later \code{PREDICTS} formatted object need to be supplied via ...
#' @param path_species A [`character`] pointing to a `PREDICTS` formatted occurrence database.
#' @param path_sites A [`character`] pointing to a `PREDICTS` formatted site database.
#' @param version A [`character`] denoting the version of the PREDICTS database to use. Has to be either
#' \code{"2016"}, \code{"2022"} or \code{"custom"}. In the latter case it is important that formatting is consistent
#' with earlier version of the PREDICTS database (Default: \code{"2016"}).
#' @param output A optional [`character`] name to store the underlying database. If set to \code{NULL} the database is only created in
#' a temporary directory (Default: \code{NULL}).
#' @param ... Other non-specified parameters.
#' @references 
#' * 1. Hudson, L.N., Newbold, T., Contu, S., Hill, S.L.L., Lysenko, I., De Palma, A., Phillips, H.R.P., Senior, R.A., Bennett, D.J., Booth, H., et al. (2014). The PREDICTS database: a global database of how local terrestrial biodiversity responds to human impacts. Ecol. Evol. 4, 4701–4735. 10.1002/ece3.1303.
#' * 2. Hudson, L.N., Newbold, T., Contu, S., Hill, S.L.L., Lysenko, I., De Palma, A., Phillips, H.R.P., Alhusseini, T.I., Bedford, F.E., Bennett, D.J., et al. (2017). The database of the PREDICTS (Projecting Responses of Ecological Diversity In Changing Terrestrial Systems) project. Ecol. Evol. 7, 145–188. 10.1002/ece3.2579.
#' @keywords database, preparation
#' @examples
#' \dontrun{
#'   path_species <- "../../100_rawdata/PREDICTS/2016/database.rds"
#'   path_sites <- "../../100_rawdata/PREDICTS/2016/sites.rds"
#'   obj <- create_hib(path_species, path_sites)
#'   obj
#' }
#' @export
create_hib <- function(path_species, path_sites, version = "2016", output = NULL, ...){
  assertthat::assert_that(
    is.character(path_species),
    file.exists(path_species),
    is.character(path_sites),
    file.exists(path_sites),
    is.character(version),
    is.null(output) || is.character(output)
  )

  # Check that version is correct
  version <- match.arg(version, c("2016", "2022", "custom"), several.ok = TRUE)
  
  # If output is null, create one
  if(is.null(output)){
    output <- paste0(tempfile(),".rds")
  }
  if(!assertthat::has_extension(output, "rds")){
    output <- paste0(output, ".rds")
  }
  
  # ------ #
  message("Create database object...")
  obj <- Database$new(
    version = version,
    paths = list(sites = path_sites,
                 species = path_species,
                 path_output = output)
  )
  
  # Save the output (using the output path specified above)
  saveRDS(obj, output)
  assertthat::assert_that(file.exists(output),
                          msg = "Writing a copy of the database failed...")
  
  return(obj)
}

#' Load a previously created database object
#' 
#' @description 
#' The function [`create_hib`] allows to save the output of the formatted database object in a specific
#' place. This function can load such a previously created database, allowing interaction with the output. 
#' @details 
#' To be run after [`create_hib`]. Expecting a file with \code{rds} ending.
#' @param path_input A [`character`] pointing to a filename.
#' @keywords database, preparation
#' @examples 
#' \dontrun{
#'   obj <- load_hib(path_input)
#' } 
#' @export
load_hib <- function(path_input){
  assertthat::assert_that(
    is.character(path_input),
    file.exists(path_input),
    assertthat::has_extension(path_input, "rds")
  )
  # Try and read in the database
  obj <- readr::read_rds(path_input)
  
  #Make some sensitivity checks that the object is correct
  obj$checkself()
  
  return(obj)
}