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
                          # dataset_name = "", ...){
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
    # is.character(dataset_name)
  )

  # TODO: replace version with dataset_name
  # if (dataset_name == ""){
  #   dataset_name <- tools::file_path_sans_ext(basename(path_species))
  # }

  # TODO: 'version' feels restrictive/wrong - can load non-PREDICTS databases
  # Maybe replace with 'dataset_name'
  # 
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
  
  # ------ #
  # ??
  # message(crayon::green("Creating HibDatabase object...\n"))
  
  # dataset_name
  names(path_species) <- version
  names(path_sites) <- version
  
  obj <- HibDatabase$new(
    version = version,
    # dataset_name = dataset_name,
    paths = list(sites = path_sites,
                 species = path_species)
  )
  
  return(obj)
}


