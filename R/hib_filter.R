#' Filter the database object by site(s)
#' 
#' @description 
#' This function returns a database object containing only sites whose characteristics 
#' are retained
#' @param obj A [`hib_database`] object.
#' @param condition An expression identifying how to filter the sites dataframe. 
#' @keywords database filtering
#' @examples
#' \dontrun{
#'   hib_filter_species(hib_obj, UN_region == "Europe")
#' }
#' @export
hib_filter_sites <- function(obj, condition){
  # str_cond <- deparse(rlang::quo_get_expr(dplyr:::dplyr_quosures({{condition}})[[1]]))
  # str_cond <- gsub("\"", "'", str_cond)
  
  assertthat::assert_that(
    inherits(obj, "Database")
    )

  # TODO: Is there a way to check the condition is a condition
  # print({{condition}})

  # print(is.character({{condition}}))
  # print(substitute(condition))
  # print(dplyr:::dplyr_quosures({{condition}})[[1]])
  # print(deparse(condition))
  # TODO: if condition is an expression, use deparse
  # cond_type <- try(inherits(condition, "expression"))
  # print(cond_type)

  str_cond <- dplyr:::dplyr_quosures({{condition}})[[1]] |>
    rlang::quo_get_expr() |>
    deparse() |>
    trimws() |>
    paste(collapse = " ")
  
  str_cond <- gsub("\"", "'", str_cond)
  str_cond <- gsub("  ", "", str_cond)
  
  message(crayon::green(paste("Filtering database using:", 
                              str_cond,
                              "...\n")))
  
  sites_upd <- dplyr::filter(obj$sites, {{ condition }} )
  species_upd <- dplyr::filter(obj$species, SSBS %in% sites_upd$SSBS)
  
  assertthat::assert_that( 
    nrow(sites_upd) > 0,
    nrow(species_upd) > 0
  )
  # TODO: any need for further assertions??
  # TODO: could compare nrow before and after - throw a warning if no reduction
  # %%%
  nrow_site_orig <- nrow(obj$sites)
  nrow_site_new <- nrow(sites_upd)

  if (nrow_site_new == nrow_site_orig){
    message(crayon::yellow(paste("Warning: the number of rows in the sites data",
                                 "has not changed following filtering.")))
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
  obj$last_update <- lubridate::today()
  
  obj$checkself()
  # Print updated summary of database
  # message(crayon::green("Filtered database summary:"))
  # print(obj)
  # message("\n")
  # Re-calculate effort-adjusted measurements? - ensure users filter first..!!
  
  return(obj)
}


#' Filter the database object by species or groups of species
#' 
#' @description 
#' This function returns a database object containing only species (and associated sites) whose characteristics 
#' are retained
#' @param obj A [`hib_database`] object.
#' @param condition An expression identifying how to filter the species dataframe
#' @keywords database, filtering
#' @examples
#' \dontrun{
#'   hib_filter_species(hib_obj, Class == "Insecta")
#' }
#' @export
hib_filter_species <- function(obj, condition){
  # e.g. by taxa: "Kingdom", "Phylum",  "Class", "Order", "Family", "Genus"
  # str_cond <- deparse(rlang::quo_get_expr(dplyr:::dplyr_quosures({{condition}})[[1]]))
  # str_cond <- gsub("\"", "'", str_cond)
  assertthat::assert_that(
    inherits(obj, "Database")
    )
  # TODO: Is there a way to check the condition is a condition

  str_cond <- dplyr:::dplyr_quosures({{condition}})[[1]] |>
    rlang::quo_get_expr() |>
    deparse() |>
    trimws() |>
    paste(collapse = " ")
  
  str_cond <- gsub("\"", "'", str_cond)
  str_cond <- gsub("  ", "", str_cond)
  
  message(crayon::green(paste("Filtering database using:", 
                              str_cond,
                              "...\n")))
  
  species_upd <- dplyr::filter(obj$species, {{ condition }} )
  
  assertthat::assert_that( nrow(species_upd) > 0 )
  # TODO: any need for further assertions??
  
  # Re-calculate site summary stats...
  sites_upd <- species_upd |>
    dplyr::group_by(SSBS) |>
    dplyr::summarise(Higher_taxa = paste0(unique(Higher_taxon), collapse = ","),
                     N_samples = dplyr::n()) |>
    dplyr::ungroup()
  
  sites_orig_sub <- dplyr::filter(obj$sites, SSBS %in% sites_upd$SSBS)
  
  # TODO: what if Higher_taxa/N_samples are not in sites, 
  sites_upd <- sites_orig_sub |>
    dplyr::select(-c(Higher_taxa, N_samples)) |>
    dplyr::left_join(sites_upd, by = "SSBS")
  
  nrow_spp_orig <- nrow(obj$species)
  nrow_spp_new <- nrow(species_upd)
  if (nrow_spp_new == nrow_spp_orig){
    message(crayon::yellow(paste("Warning: the number of rows in the species data",
                                 "has not changed following filtering.")))
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
  obj$last_update <- lubridate::today()
  
  obj$checkself()
  
  # message(crayon::green("Filtered database summary:"))
  # print(obj)
  # message("\n")
  # Re-calculate effort-adjusted measurements? - ensure users filter first..!!
  
  return(obj)
}


# TODO
# What about filtering sites if any spp meets a criteria...
# Is this sensible?


