#' Create a database object from input files
#' 
#' @description 
#' This function creates an [`R6`] `HibDatabase` object containing the loaded  
#' \code{PREDICTS}-like database.
#' It assumes that a species database file is present, with a site database file 
#' being optional. 
#' All files are loaded and compiled into a single database object that can be 
#' stored for future use (see [hib_save()]).
#' 
#' This object is then used in other functions of this package to prepare or filter the
#' database and enable ecological impact analyses. 
#' @note 
#' Other or later \code{PREDICTS} formatted object need to be supplied via ... ??
#' @param path_species A character pointing to a `PREDICTS` formatted abundance/
#' occurrence database.
#' @param path_sites Optional character argument, pointing to a `PREDICTS` 
#' formatted site database. Default value is `""`.
#' @param dataset_name A [`character`] indicating the dataset used to create the
#' HibDatabase. Default is `""`. 
#' If `""` is supplied, the parameter will be internally updated to the filename 
#' associated with \code{path_species}.
#' @param ... Other non-specified parameters.
#' @returns An object of class `HibDatabase`, containing site and species information 
#' which can be further modified for subsequent ecological analysis. 
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
hib_create_db <- function(path_species, 
                          path_sites = "", 
                          dataset_name = "", 
                          ...){

  # Assertions
  assertthat::assert_that(
    is.character(path_species),
    file.exists(path_species),
    is.character(path_sites),
    (file.exists(path_sites) || path_sites == ""),
    is.character(dataset_name)
  )

  # dataset_name
  if (dataset_name == ""){
    dataset_name <- tools::file_path_sans_ext(basename(path_species))
  }
 
  names(path_species) <- dataset_name
  names(path_sites) <- dataset_name
  
  obj <- HibDatabase$new(
    dataset_name = dataset_name,
    paths = list(sites = path_sites,
                 species = path_species)
  )
  
  return(obj)
}


