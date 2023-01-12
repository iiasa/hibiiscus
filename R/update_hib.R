#' This function updates an existing Database object with further data
#' 
#' @description 
#' The function [`create_hib`] allows to save the output of the formatted database object in a specific
#' place. This function can update an existing database, for example by adding the 2022 or a custom PREDICTs version 
#' to it.
#' @details 
#' To be run after [`create_hib`]. Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @inheritParams create_hib
#' @keywords database, preparation
#' @examples 
#' \dontrun{
#'   obj <- load_hib(path_input)
#' } 
#' @export
update_hib <- function(obj, path_species, path_sites, version = "2016", output = NULL, ...){
    assertthat::assert_that(
      inherits(obj, "Database"),
      is.character(path_species),
      file.exists(path_species),
      is.character(path_sites),
      file.exists(path_sites),
      is.character(version),
      is.null(output) || is.character(output)
    )

    # Check that version is correct
    version <- match.arg(version, c("2016", "2022", "custom"), several.ok = TRUE)
    # Ensure that the version isn't already in the database
    assertthat::assert_that(!version %in% obj$version,
                            msg = "This version of PREDICTS is already in the database object!")
    
    # To checks of the existing database
    obj$checkself()
    
    # Merge the new data with the existing database and update the object!
    # TODO: Also update version number and last_updated to it
    stop("To be coded!")
        
    return(obj)
}