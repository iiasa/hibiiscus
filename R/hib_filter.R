#' hib_filter_sites: Filter the database object by site-linked variables
#' 
#' @description 
#' This function returns a database object containing only sites (and associated
#' species) whose characteristics are retained.
#' 
#' @param obj A [`HibDatabase`] object.
#' @param condition An expression identifying how to filter the sites data.frame, 
#' as would be passed to `dplyr::filter`.
#' 
#' @return `HibDatabase`
#' 
#' @keywords database, filtering
#' 
#' @examples
#' \dontrun{
#'   hib_obj <- hib_filter_sites(hib_obj, UN_region == "Europe")
#' }
#' @export
hib_filter_sites <- function(obj, condition){
  
  assertthat::assert_that(
    inherits(obj, "HibDatabase")
    )

  # TODO: Is there a way to check the condition is a condition
  # TODO: if condition is an expression, use deparse
  
  str_cond <- dplyr:::dplyr_quosures({{condition}})[[1]] |>
    rlang::quo_get_expr() |>
    deparse() |>
    trimws() |>
    paste(collapse = " ")
  
  str_cond <- gsub("\"", "'", str_cond)
  str_cond <- gsub("  ", "", str_cond)
  
  message(crayon::green(paste("\nFiltering HibDatabase using:", 
                              str_cond,
                              "...\n")))
  
  sites_upd <- dplyr::filter(obj$sites, {{ condition }} )
  species_upd <- dplyr::filter(obj$species, SSBS %in% sites_upd$SSBS)
  
  assertthat::assert_that( 
    nrow(sites_upd) > 0,
    nrow(species_upd) > 0
  )
  
  # Compare nrow before and after - throw a warning if no reduction
  # %%%
  nrow_site_orig <- nrow(obj$sites)
  nrow_site_new <- nrow(sites_upd)

  if (nrow_site_new == nrow_site_orig){
    message(crayon::yellow(paste("Warning: the number of rows in the sites data",
                                 "has not changed following filtering.\n")))
  }
  # %%%

  # Update obj
  obj$species <- species_upd
  obj$sites <- sites_upd
  # Append filtering flag to database object
  filter_idx <- length(obj$processing_flags$filtering_steps) + 1
  
  obj$processing_flags$filtering_steps[[filter_idx]] <- paste0(as.character(filter_idx),
                                                               ": ",
                                                               str_cond)
  obj$last_update <- lubridate::now()
  
  obj$checkself()
  
  return(obj)
}


#' hib_filter_species: Filter the database object by species-linked features.
#' 
#' @description 
#' This function returns a database object containing only species (and associated
#' sites) whose characteristics are retained.
#' 
#' @param obj A [`HibDatabase`] object.
#' @param condition An expression identifying how to filter the species data.frame,
#' as would be passed to `dplyr::filter`.
#' 
#' @return `HibDatabase`
#' 
#' @keywords database, filtering
#' @examples
#' \dontrun{
#'   hib_obj <- hib_filter_species(hib_obj, Class == "Insecta")
#' }
#' @export
hib_filter_species <- function(obj, condition){

  assertthat::assert_that(
    inherits(obj, "HibDatabase")
    )
  # TODO: Is there a way to check the condition is a condition

  str_cond <- dplyr:::dplyr_quosures({{condition}})[[1]] |>
    rlang::quo_get_expr() |>
    deparse() |>
    trimws() |>
    paste(collapse = " ")
  
  str_cond <- gsub("\"", "'", str_cond)
  str_cond <- gsub("  ", "", str_cond)
  
  message(crayon::green(paste("\nFiltering HibDatabase using:", 
                              str_cond,
                              "...\n")))
  
  species_upd <- dplyr::filter(obj$species, {{ condition }} )
  
  assertthat::assert_that( nrow(species_upd) > 0 )
  
  # Use species_upd to filter sites df
  sites_upd <- species_upd |>
    dplyr::group_by(SSBS) |>
    dplyr::summarise(Higher_taxa = paste0(unique(Higher_taxon), collapse = ","),
                     N_samples = dplyr::n()) |>
    dplyr::ungroup()
  
  sites_orig_sub <- dplyr::filter(obj$sites, SSBS %in% sites_upd$SSBS)
  
  sites_upd <- sites_orig_sub |>
    dplyr::select(-c(Higher_taxa, N_samples)) |>
    dplyr::left_join(sites_upd, by = "SSBS")
  
  nrow_spp_orig <- nrow(obj$species)
  nrow_spp_new <- nrow(species_upd)
  if (nrow_spp_new == nrow_spp_orig){
    message(crayon::yellow(paste("Warning: the number of rows in the species data",
                                 "has not changed following filtering.\n")))
  }
  # %%%

  # Update obj
  obj$species <- species_upd
  obj$sites <- sites_upd
  # Append filtering flag to database object
  filter_idx <- length(obj$processing_flags$filtering_steps) + 1
  
  obj$processing_flags$filtering_steps[[filter_idx]] <- paste0(as.character(filter_idx),
                                                               ": ",
                                                               str_cond)
  obj$last_update <- lubridate::now()
  
  obj$checkself()
  
  return(obj)
}


# TODO
# Study-level filters?...
# What about filtering sites if any spp meets a criteria...


