#' Filter the database object by site(s)
#' 
#' @description 
#' This function returns a database object containing only sites whose characteristics 
#' are retained
#' @param obj A [`hib_database`] object.
#' @param condition A [`conditional`] statement identifying how to filter the sites dataframe
#' @keywords database, filtering
#' #' @examples
#' \dontrun{
#'   hib_filter_species(hib_obj, UN_region == "Europe")
#' }
hib_filter_sites <- function(obj, condition){
  str_cond <- deparse(rlang:::quo_get_expr(dplyr:::dplyr_quosures({{condition}})[[1]]))
  str_cond <- gsub("\"", "'", str_cond)
  
  message(crayon::green(paste("Filtering database using:", 
                              str_cond,
                              "...\n")))
  
  sites_upd <- dplyr::filter(obj$sites, {{ condition }} )
  species_upd <- dplyr::filter(obj$species, SSBS %in% sites_upd$SSBS)
  
  assertthat::assert_that( 
    nrow(sites_upd) > 0,
    nrow(species_upd) > 0)
  # TODO: any need for further assertions??
  
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
  message(crayon::green("Filtered database summary:"))
  print(obj)
  message("\n")
  # Re-calculate effort-adjusted measurements? - ensure users filter first..!!
  
  return(obj)
}


#' Filter the database object by species or groups of species
#' 
#' @description 
#' This function returns a database object containing only species (and associated sites) whose characteristics 
#' are retained
#' @param obj A [`hib_database`] object.
#' @param condition A [`conditional`] statement identifying how to filter the species dataframe
#' @keywords database, filtering
#' @examples
#' \dontrun{
#'   hib_filter_species(hib_obj, Class == "Insecta")
#' }
hib_filter_species <- function(obj, condition){
  # e.g. by taxa: "Kingdom", "Phylum",  "Class", "Order", "Family", "Genus"
  str_cond <- deparse(rlang:::quo_get_expr(dplyr:::dplyr_quosures({{condition}})[[1]]))
  str_cond <- gsub("\"", "'", str_cond)
  
  message(crayon::green(paste("Filtering database using:", 
                              str_cond,
                              "...\n")))
  
  species_upd <- dplyr::filter(obj$species, {{ condition }} )
  
  assertthat::assert_that( nrow(species_upd) > 0 )
  # TODO: any need for further assertions??
  
  # Re-calculate site summary stats...
  sites_upd <- species_upd %>%
    dplyr::group_by(SSBS) %>%
    dplyr::summarise(Higher_taxa = paste0(unique(Higher_taxon), collapse = ","),
                     N_samples = n()) %>%
    dplyr::ungroup()
  
  sites_orig_sub <- dplyr::filter(obj$sites, SSBS %in% sites_upd$SSBS)
  
  sites_upd <- sites_orig_sub %>%
    dplyr::select(-c(Higher_taxa, N_samples)) %>%
    dplyr::left_join(sites_upd, by = "SSBS")
  
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
  
  message(crayon::green("Filtered database summary:"))
  print(obj)
  message("\n")
  # Re-calculate effort-adjusted measurements? - ensure users filter first..!!
  
  return(obj)
}

# subset(db_2022$species, Class == "Mammalia" & Realm == "Australasia")
# test <- hib_filter_species_alt(db_2022, Class == 'Mammalia' & Realm == 'Australasia')
# print(test)
# deparse(text = test)
# 
# deparse(test)


# hib_filter_species <- function(obj, column_nam, vals_to_keep){
# stop("TBD")
# Re-generate site summary data - e.g. N samples may have changed...spp considered
# Re-calculate effort-adjusted measurements? - ensure users filter first..!!
# Append filtering flag to database object
# taxa-focussed filter: ....
# }


# Basically need to re-write a daplyr::filter function
# the after filtering species/sites, filters the other too
# and then recalcs site metics/effort_corrected_measurements...
# need a style of hib_filter_species(obj, filter_formula)
# where filter_formula is something like Family == "X", or Family != "X"


# What about filtering sites if any spp meets a criteria...
# Is this sensible?


