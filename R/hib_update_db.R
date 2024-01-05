#' This function updates an existing Database object with further data
#' 
#' @description 
#' This function can update an existing database, for example by adding the 2022 
#' or a custom PREDICTS version to it.
#' @details 
#' Expects a [`HibDatabase`] object.
#' @param obj A [`HibDatabase`] object.
#' @inheritParams hib_create_db
#' @keywords database, preparation
#' @examples 
#' \dontrun{
#'   obj <- update_hib(obj1, obj2)
#' } 
#' @export
hib_update_db <- function(obj, path_species, path_sites = "", 
                          # dataset_name = "", 
                       version = "2022", custom_version_name = "", ...){
# update_hib

  assertthat::assert_that(
    inherits(obj, "HibDatabase"),
    is.character(path_species),
    file.exists(path_species),
    (is.character(path_sites) & file.exists(path_sites) ) || (is.character(path_sites) & path_sites == ""),
    is.character(version),
    is.character(custom_version_name)
    # is.character(dataset_name)
    
  )
  # TODO: replace version with dataset_name
  # if (dataset_name == ""){
  #   dataset_name <- tools::file_path_sans_ext(basename(path_species))
  # }

    # Check that version is correct
    version <- match.arg(version, c("2016", "2022", "custom"), several.ok = F)
    
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
    
    # Ensure that the version isn't already in the database
    assertthat::assert_that(!version %in% obj$version,
                            msg = "This version of PREDICTS is already in the HibDatabase object!")
    # This dataset (dataset_name) is ...


    # Checks of the existing database
    obj$checkself()
    
    message(crayon::green("Updating HibDatabase ...\n"))
    
    # dataset_name
    names(path_species) <- version
    names(path_sites) <- version
    # Create new object
    new_obj <- HibDatabase$new(
      version = version,
      # dataset_name = dataset_name,
      paths = list(sites = path_sites,
                   species = path_species)#,
      # path_output = output)
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
    
    # Update version info and update info
    # dataset_name
    old_version <- obj$version

    obj$version <- paste(obj$version, version, sep = "; ")
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
    obj$processing_flags$database_updates[[join_idx]] <- paste0(as.character(join_idx),
                                                                  ": ",
                                                                  "Updatated database - Added '",
                                                                  version,
                                                                  "' to '",
                                                                  old_version, 
                                                                  "'.")
    # %%%%%

    # Checks of the updated database
    obj$checkself()

    return(obj)
}