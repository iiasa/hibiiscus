#' hib_create_lu_legend: Create a 'new' land-use legend.
#' 
#' @description 
#' This function creates a [`R6`] [`HibLULegend`] object containing the provided 
#' land-use legend (i.e., linking sites, SSBS, to land-uses/intensities). 
#' 
#' @details
#' The provided data.frame must include the following columns:
#' `SSBS`, `Predominant_land_use`, `Use_intensity` and `LUI`.
#' 
#' @param input Either a `character` pointing to a saved rds file, or a 
#' `data.frame`-like object within the environment.
#' 
#' @return `HibLULegend`
#' 
#' @keywords database, land-use, legend
#' 
#' @seealso [`HibLULegend`]
#' 
#' @examples
#' \dontrun{
#'   lu_leg <- hib_create_lu_legend("../Data/landuse_dataframe.rds")
#' }
#' @export
hib_create_lu_legend <- function(input){
  
  assertthat::assert_that(
    is.character(input) || is.data.frame(input)
  )
  if (is.character(input)){
    assertthat::assert_that(
      file.exists(input),
      assertthat::has_extension(input, "rds"))
  }
  
  leg <- HibLULegend$new(input = input,
                         input_chr = deparse(substitute(input)) )
  
  return(leg)
}



#' hib_set_lu_legend: Set a `HibLULegend` as the LU legend in a `HibDatabase`
#'
#' @description
#' This function adds the `HibLULegend` to the `HibDatabase`, modifying both the
#' `sites` and `species` data.frames.
#' 
#' @details
#' This function updates (and therefore overwrites) the previous 
#' `Predominant_land_use`, `Use_intensity` and `LUI` columns in both the `sites` 
#' and `species` data.frames.
#'
#' @param obj A `HibDatabase`.
#' @param leg A `HibLULegend`.
#' 
#' @return `HibDatabase`
#' 
#' @seealso 
#' [HibDatabase]
#' [HibLULegend]
#' 
#' @keywords database, land-use, legend
#' 
#' @examples
#' \dontrun{
#'   obj <- hib_set_lu_legend(obj, leg)
#' }
#' @export
hib_set_lu_legend <- function(obj, leg){
  
  assertthat::assert_that(
    inherits(obj, "HibDatabase"),
    inherits(leg, "HibLULegend") 
  )
  
  message(crayon::green("\nSetting land-use legend for HibDatabase...\n"))
  
  # Make some sensitivity checks that the inputs are okay
  obj$checkself()
  leg$checkself()
  
  # check for column names already in sites (Predominant_land_use, Use_intensity, LUI)
  cols_to_drop_sites <- colnames(obj$sites)[colnames(obj$sites) %in% c("Predominant_land_use", "Use_intensity", "LUI")]
  # Drop cols from sites
  obj$sites <- obj$sites |>
    dplyr::select(!dplyr::all_of(cols_to_drop_sites)) |>
    # left join new legend info
    dplyr::left_join(leg$lu_df)
  
  cols_to_drop_spp <- colnames(obj$species)[colnames(obj$species) %in% c("Predominant_land_use", "Use_intensity", "LUI")]
  # Drop cols from species
  obj$species <- obj$species |>
    dplyr::select(!dplyr::all_of(cols_to_drop_spp)) |>
    # left join new legend info
    dplyr::left_join(leg$lu_df)

  # Add in a legend flag
  obj$processing_flags$lu_legend_source <- leg$input_data
  
  obj$last_update <- Sys.time()
  
  return(obj)
  
}

