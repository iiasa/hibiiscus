#' hib_update_db: Update an existing \code{HibDatabase} with further data
#' 
#' @description 
#' This function updates an existing \code{HibDatabase}, for example by adding 
#' the 2022 version of the PREDICTS data to the 2016 version.
#' 
#' @details 
#' Newly provided data (from files) will be appended to the previously created
#' \code{HibDatabase}. Species data will be added to the species data.frame, and 
#' site data to the site data.frame.
#' Meta-data/descriptors in the \code{HibDatabase} will be updated to reflect 
#' the additions.
#' Note: Prior filtering/processing of the initial \code{HibDatabase} will not 
#' be applied to the newly provided data. Please load all data of interest into 
#' the \code{HibDatabase} prior to further processing.
#' 
#' @param obj A [`HibDatabase`] object.
#' @inheritParams hib_create_db
#' @returns The updated [`HibDatabase`] object, with the supplied dataset added to its site and species data.
#' @inherit hib_create_db references
#' 
#' @keywords database, preparation
#' @examples 
#' \dontrun{
#'   path_species <- "rawdata/PREDICTS/2016/database.rds"
#'   path_sites <- "rawdata/PREDICTS/2016/sites.rds"
#'   hibDB <- hib_create_db(path_species, path_sites, "predicts_2016")
#'   path_species_upd <- "rawdata/PREDICTS/2022/database2ndrelease.rds"
#'   hibDB <- update_hib(hibDB, path_species_upd, "", "predicts_2022")
#'   hibDB
#' } 
#' @export
hib_update_db <- function(obj, 
                          path_species, 
                          path_sites = "", 
                          dataset_name = ""){

  assertthat::assert_that(
    inherits(obj, "HibDatabase"),
    is.character(path_species),
    file.exists(path_species),
    (is.character(path_sites) & file.exists(path_sites) ) || (is.character(path_sites) & path_sites == ""),
    is.character(dataset_name)
    
  )
  
  if (dataset_name == ""){
    dataset_name <- tools::file_path_sans_ext(basename(path_species))
  }
  
  # Ensure that the version isn't already in the database
  assertthat::assert_that(!dataset_name %in% obj$dataset_name,
                          msg = "This dataset (dataset_name) is already in the HibDatabase object!")
  
  # Checks of the existing database
  obj$checkself()
  
  message(crayon::green("Updating HibDatabase ..."))
  
  # dataset_name
  names(path_species) <- dataset_name
  names(path_sites) <- dataset_name
  
  # Create new object
  new_obj <- HibDatabase$new(
    dataset_name = dataset_name,
    paths = list(sites = path_sites,
                 species = path_species)
  )
  
  # Assert that no sites are already in current version...
  # If they are, throw a warning and drop them
  new_sites_to_rm <- unique(new_obj$sites$SSBS)[unique(new_obj$sites$SSBS) %in% unique(obj$sites$SSBS)]
  if (length(new_sites_to_rm) > 0 ){
    message(crayon::yellow(paste0("Warning: ", length(new_sites_to_rm),
                                  " sites in the new data are alrady present in the HibDatabase!",
                                  "\n",
                                  "These sites will not be updated.\n")))
    print(new_sites_to_rm)
    
    new_obj$species <- new_obj$species |>
      dplyr::filter(!(SSBS %in% new_sites_to_rm))
    new_obj$sites <- new_obj$sites |>
      dplyr::filter(!(SSBS %in% new_sites_to_rm))
  }
  
  # Add species/sites to original object
  obj$species <- dplyr::bind_rows(obj$species,
                                  new_obj$species)    
  
  obj$sites <- dplyr::bind_rows(obj$sites,
                                new_obj$sites)    
  
  # Update dataset_name info and update info
  old_dat_nm <- obj$dataset_name
  
  obj$dataset_name <- paste(obj$dataset_name, dataset_name, sep = "; ")
  obj$last_update <- lubridate::now()
  
  # Add paths...
  obj$paths$species <- c(obj$paths$species, new_obj$paths$species)
  obj$paths$sites <- c(obj$paths$sites, new_obj$paths$sites)
  
  # TODO
  # Potential extra issues - previous filtering of old version/div object?
  # Can rm diversity?, but filtering etc can't be undone
  
  # Add in a processing flag for update:
  # "Updated HibDatabase: added <X> data to <Y> data."
  # %%%%%
  
  # Append update flag to database object
  upd_idx <- length(obj$processing_flags$database_updates) + 1
  
  # dataset_name
  obj$processing_flags$database_updates[[upd_idx]] <- paste0(as.character(upd_idx),
                                                             ": ",
                                                             "Updatated database - Added '",
                                                             dataset_name,
                                                             "' to '",
                                                             old_dat_nm, 
                                                             "'.")
  # %%%%%
  
  # Checks of the updated database
  obj$checkself()
  
  return(obj)
}