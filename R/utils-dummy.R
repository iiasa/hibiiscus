#' Internal function to create dummy data for PREDICTS style models
#' 
#' @description
#' The package itself should not contain or reference any external datasets or 
#' or their paths (which are system specific anyway).
#' Rather a dummy dataset with the same format should be used, also allowing 
#' us to create tests for the package that are not reliant on external data.
#' 
#' @details
#' This function creates dummy data in the format of the PREDICTS database. 
#' A separate function (TBD) can be used to create external information for 
#' projection or integration
#' 
#' @param what [`character`] stating what exactly to create. Available options are
#' \code{"sites"} to emulate a PREDICTS site database and \code{"species"} for a species database.
#' @param n_studies [`numeric`] on the number of studies to be created.
#' @param n_sites [`numeric`] Number of sites per study to be created.
#' @param lulegend [`character`] of the land use legend to be generated for the dummy data.
#' @examples
#' \dontrun{
#'   out <- create_dummy_data()
#'   # Save in hibiscus data folder (lazy loading enabled)
#'   save(out, "data/dummydata.rda")
#' }
#' @import stringi stri_rand_strings
#' @keywords internal
#' @author Martin Jung
#' @noRd
create_dummy_data <- function(what = "sites",
                              n_studies = 20, n_sites = 10,
                              lulegend = c("Primary", "Secondary", "Cropland", "Urban")
                              ){
  # Check input parameters
  assertthat::assert_that(
    is.character(what),
    is.numeric(n_studies), n_studies >= 3,
    is.numeric(n_sites), n_sites >= 2,
    is.character(lulegend), length(lulegend) >= 2
  )
  # Match to available options
  what <- match.arg(what, c("sites", "species"),several.ok = FALSE)
  
  # Make the land-use legend an factor, with the first level being the reference
  lulegend <- factor(lulegend, levels = c(lulegend[1], lulegend[2:length(lu_legend)]))
  
  # Normalize function
  norm <- function(y) (y - min(y)) / (max(y) - min(y)) 
  
  # Create a new data.frame container
  if(what == "sites"){
    
    # Create ids of studies
    ns <- stringi::stri_rand_strings(n = n_studies,length = 7)
    assertthat::assert_that(length(unique(ns)) == n_studies)
    lui <- c("Minimal use", "Light use", "Intense use")
   
    out <- tidyr::crossing(
      SS = paste0(ns, " ", 1),
      SSBS = 1:n_sites
    ) |> 
      dplyr::mutate(SSBS = paste0(SS, " ", SSBS),
                    # TODO: Add blocking scheme
                    Source_id = substr(SS, start = 0, stop = 7),
                    Predominant_land_use = sample(lulegend, size = n_sites * n_studies, replace = TRUE),
                    Use_intensity = sample(lui, size = n_sites * n_studies, replace = TRUE, prob = c(1,.5,.2))
                    )

    # Other study settings
    # Normalize and invert
    out$Rescaled_sampling_effort <- 1- norm(rexp(nrow(out),rate = 2))
    
    # Generate dummy coordinates
    # TODO:
    
    ## Generate some random measurements ##
    out$Species_richness <- rpois(n = nrow(out), 5) 
    out$LogAbundance <- rnorm(n = nrow(out),10, .5)
      
  } else {
    stop("Not yet implemented!")
  }
  return(out)
}
