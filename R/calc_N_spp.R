#' calc_N_spp
#' 
#' @description 
#' Calculate number of species contributing to a model (internal function).
#' 
#' @details 
#' This function calculates the number of species contributing to a model fit 
#' based on the metric (model response) of interest. Intended for internal use.
#' Calls \code{get_n_spp_taxa}.
#' 
#' @param obj A HibDatabase object.
#' @param mod_df A data.frame object created as part of data preparation for
#' model fitting.
#' @param mod_resp A [`character`] object indicating the modelled response variable. - NOT USED
#' @return list
#' 
#' @keywords database, model, processing
#' 
#' @seealso
#' `get_n_spp_taxa()`
#' 
#' @examples 
#' \dontrun{
#'   n_spp <- calc_N_spp(obj, mod_df)
#' } 
#' 
#' @noRd
calc_N_spp <- function(obj,
                       mod_df,
                       mod_resp){
  
  assertthat::assert_that(
    inherits(obj, "HibDatabase"),
    is.data.frame(mod_df)
    )
  
  # Check if measurement or effort corrected 
  # When would one be >0 and the other not
  if (obj$processing_flags$measures_effort_corrected == T){
    measure_var <- "Effort_corrected_measurement"
  } 
  if (obj$processing_flags$measures_effort_corrected == F){
    measure_var <- "Measurement"
  }

  if (all(c("SSBS1", "SSBS2") %in% colnames(mod_df))){
    ssbs <- unique(c(mod_df$SSBS1, mod_df$SSBS2))
  } else {
    ssbs <- unique(mod_df$SSBS)
  }

  # Filter obj$species to SSBS/SSBS1/SSBS2 in mod_df
  # then calc N_spp/N_taxa
  df <- obj$species |>
    dplyr::filter(SSBS %in% ssbs)
  out <- get_n_spp_taxa(df, measure_var)

  return(out)
  
}

