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
#' @param output A optional [`character`] name to store the underlying database (Default: \code{NULL}).
#' @param ... Other non-specified parameters.
#' @references 
#' * 1. Hudson, L.N., Newbold, T., Contu, S., Hill, S.L.L., Lysenko, I., De Palma, A., Phillips, H.R.P., Senior, R.A., Bennett, D.J., Booth, H., et al. (2014). The PREDICTS database: a global database of how local terrestrial biodiversity responds to human impacts. Ecol. Evol. 4, 4701–4735. 10.1002/ece3.1303.
#' * 2. Hudson, L.N., Newbold, T., Contu, S., Hill, S.L.L., Lysenko, I., De Palma, A., Phillips, H.R.P., Alhusseini, T.I., Bedford, F.E., Bennett, D.J., et al. (2017). The database of the PREDICTS (Projecting Responses of Ecological Diversity In Changing Terrestrial Systems) project. Ecol. Evol. 7, 145–188. 10.1002/ece3.2579.
#' @keywords database, preparation
#' @examples
#' \dontrun{
#'   path_species <- "../../../PhD/Data/PREDICTS_v1/database.rds"
#'   path_sites <- "../../../PhD/Data/PREDICTS_v1/sites.rds"
#'   obj <- create_hib(path_species, path_sites)
#'   obj
#' }
#' @export
create_hib <- function(path_species, path_sites, output = NULL, ...){
  assertthat::assert_that(
    is.character(path_species),
    file.exists(path_species),
    is.character(path_sites),
    file.exists(path_sites),
    is.null(output) || is.character(output)
  )
  
  # If output is null, create one
  if(is.null(output)){
    output <- paste0(tempfile(),".rds")
  }
  if(!assertthat::has_extension(output, "rds")){
    output <- paste0(output, ".rds")
  }
  
  # ------ #
  message("Reading input datasets...")
  # Now load both datasets
  sites <- readr::read_rds(path_sites)
  species <- readr::read_rds(path_species)
  # ------ #
  # Security checks
  assertthat::assert_that(
    assertthat::has_name(sites, "Source_ID"),
    assertthat::has_name(sites, "SS"),
    assertthat::has_name(sites, "Diversity_metric"),
    assertthat::has_name(sites, "SSBS"),
    assertthat::has_name(sites, "Predominant_land_use"),
    # And the species data
    assertthat::has_name(species, "Source_ID"),
    assertthat::has_name(species, "SS"),
    assertthat::has_name(species, "SSBS"),
    assertthat::has_name(species, "Best_guess_binomial"),
    assertthat::has_name(species, "Measurement")
  )
  # Finally, check that all sites are in species
  assertthat::assert_that(
    all( unique(sites$Source_ID) %in% unique(species$Source_ID) )
  )
  # ------ #
  message("Create database object...")
  obj <- Database$new(
    sites = sites,
    species = species
    )
  
  # Save the output
  saveRDS(obj, output)
  
  return(obj)
}

#' Load a previously created database object
#' 
#' @description 
#' The function [`create_hib`] allows to save the output of the formatted database object in a specific
#' place. This function can load previously created files. 
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
  obj <- readr::read_rds(path_input)
  return(obj)
}