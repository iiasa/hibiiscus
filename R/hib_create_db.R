#' Create a database object
#' 
#' @description 
#' This function creates a [`R6`] `HibDatabase` object containing the loaded  
#' \code{PREDICTS} database.
#' It assumes that a species database is present, with a site database being optional. 
#' All files are loaded and compiled in a single database that can be stored for future use.
#' 
#' This object is then used in other functions of this package to prepare or filter the
#' database and enable ecological impact analyses. 
#' @note 
#' Other or later \code{PREDICTS} formatted object need to be supplied via ...
#' @param path_species A [`character`] pointing to a `PREDICTS` formatted occurrence database.
#' @param path_sites A [`character`] pointing to a `PREDICTS` formatted site database, or "" otherwise.
#' @param version A [`character`] denoting the version of the PREDICTS database to use. Has to be either
#' \code{"2016"}, \code{"2022"} or \code{"custom"}. In the latter case it is important that formatting is consistent
#' with earlier versions of the PREDICTS database (Default: \code{"2016"}).
#' @param custom_version_name Optional [`character`] that will help indicate the source 
#' of a custom database - the create_hib function will internally generate
#' paste("Custom - ", custom_version_name). Only used if version == "custom"
#' @param ... Other non-specified parameters.
#' @references 
#' * 1. Hudson, L.N., Newbold, T., Contu, S., Hill, S.L.L., Lysenko, I., De Palma, A., Phillips, H.R.P., Senior, R.A., Bennett, D.J., Booth, H., et al. (2014). The PREDICTS database: a global database of how local terrestrial biodiversity responds to human impacts. Ecol. Evol. 4, 4701–4735. 10.1002/ece3.1303.
#' * 2. Hudson, L.N., Newbold, T., Contu, S., Hill, S.L.L., Lysenko, I., De Palma, A., Phillips, H.R.P., Alhusseini, T.I., Bedford, F.E., Bennett, D.J., et al. (2017). The database of the PREDICTS (Projecting Responses of Ecological Diversity In Changing Terrestrial Systems) project. Ecol. Evol. 7, 145–188. 10.1002/ece3.2579.
#' @keywords database, preparation
#' @examples
#' \dontrun{
#'   path_species <- "rawdata/PREDICTS/2016/database.rds"
#'   path_sites <- "rawdata/PREDICTS/2016/sites.rds"
#'   obj <- create_hib(path_species, path_sites)
#'   obj
#' }
#' @export
hib_create_db <- function(path_species, path_sites = "", 
                       version = "2016", 
                       custom_version_name = "", ...){
# create_hib  

  # NOTE RC - dropping output path - not really needed - use save_hib instead.
  # output = NULL,
  #  #' @param output A optional [`character`] name to store the underlying database. 
  # If set to \code{NULL} the database is only created in a temporary directory (Default: \code{NULL}).
  
  # 2022 release doesn't have a site rds... - could make one? 
  # Or allow path_sites to be ""
  assertthat::assert_that(
    is.character(path_species),
    file.exists(path_species),
    is.character(path_sites),
    (file.exists(path_sites) || path_sites == ""),
    is.character(version),
    is.character(custom_version_name)
    # is.null(output) || is.character(output)
  )

  # Check that version is correct
  version <- match.arg(version, c("2016", "2022", "custom"), several.ok = F)
  
  # TODO
  if (version %in% c("2016", "2022")){
    version <- paste("PREDICTS -", version)
  }

  # if version == "custom", add a bit more detail
  if (version == "custom"){
    # either from custom_version_name
    if (custom_version_name != ""){
      version <- paste("Custom -", custom_version_name)
    }
    # or basename(paths$species), but also get rid of extension
    if (custom_version_name == ""){
      custom_version_name <- tools::file_path_sans_ext(basename(path_species))
      version <- paste("Custom -", custom_version_name)
    }
  }
  
  # # If output is null, create one
  # if(is.null(output)){
  #   output <- paste0(tempfile(),".rds")
  # }
  # if(!assertthat::has_extension(output, "rds")){
  #   output <- paste0(output, ".rds")
  # }
  
  # ------ #
  # ??
  # message(crayon::green("Creating HibDatabase object...\n"))
  
  names(path_species) <- version
  names(path_sites) <- version
  
  obj <- HibDatabase$new(
    version = version,
    paths = list(sites = path_sites,
                 species = path_species)
  )
  
  return(obj)
}


# # TODO: delete/rm below functions
# #' Save a previously created database object
# #' 
# #' @description 
# #' The function [`save_hib`] saves the output of the formatted [`HibDatabase`] object in a specific
# #' place.
# #' @details 
# #' Expecting a file with \code{rds} ending.
# #' @param obj A [`HibDatabase`].
# #' @param output_path A [`character`] pointing to a filename.
# #' @keywords save, database
# #' @export
# save_hib <- function(obj, output_path) obj$save_hib(output_path)


# #' Load a previously created database object
# #' 
# #' @description 
# #' The function [`load_hib`] loads a previously created and saved [`HibDatabase`]. 
# #' @details 
# #' Expecting a file with \code{rds} ending.
# #' @param path_input A [`character`] pointing to a filename.
# #' @keywords database, preparation
# #' @examples 
# #' \dontrun{
# #'   obj <- load_hib(path_input)
# #' } 
# #' @export
# load_hib <- function(path_input){

#   # hib_db_load ?? hib_load_db
#   assertthat::assert_that(
#     is.character(path_input),
#     file.exists(path_input),
#     assertthat::has_extension(path_input, "rds")
#   )
#   # Try and read in the database
#   message(crayon::green(paste0("Loading HibDatabase from: '", path_input, "...")))
  
#   obj <- readr::read_rds(path_input)
  
#   #Make some sensitivity checks that the object is correct
#   # TODO: check its a Database
#   # TODO: what about other hib datatypes... - models, cv, proj ??
#   assertthat::assert_that(
#     inherits(obj, "HibDatabase")
#     )
#   obj$checkself()
  
#   return(obj)
# }

# TODO: generic hib_save and hib_load function... - see hib_saveload

