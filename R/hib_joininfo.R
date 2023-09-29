#' Joins external species information to a database object
#' 
#' @description 
#' This function joins external information to a [`Database`] object. Such 
#' information could be for example prepared data on species traits.
#' @param obj A [`Database`] object.
#' @param spp_info A [`dataframe`] or similar of categorical or numeric information 
#' to be added to the database.
#' @keywords database, external
#' @export
hib_join_species <- function(obj, spp_info){

  # Assertions
  assertthat::assert_that(
    inherits(obj, "Database"),
    is.data.frame(spp_info)
    )
  # Tests
  # columns in spp_info 
  # - at least 1 must be in obj$species
  # - at least 1 must not be in obj$species
  obj_spp_cols <- colnames(obj$species)
  info_cols <- colnames(spp_info)
  
  n_matching_cols <- sum(info_cols %in% obj_spp_cols)
  
  assertthat::assert_that(
    n_matching_cols > 0 & n_matching_cols < length(info_cols),
    msg = "Columns in spp_info are not compatible with the current species data."
  )
  
  # TODO: joining data message
  message(crayon::green("Joining species-related data into hib database...\n"))
  
  obj$species <- dplyr::left_join(obj$species,
                           spp_info)  
  
  obj$last_update <- lubridate::today()
  
  # %%%
  # TODO: would be useful to flag addition of data - source name/an id..
  spp_info_chr <- deparse(substitute(spp_info))
  new_cols_chr <- info_cols[!info_cols %in% obj_spp_cols]

  # Append filtering flag to database object
  join_idx <- length(obj$processing_flags$joined_information) + 1
  
  obj$processing_flags$joined_information[[join_idx]] <- paste0(as.character(join_idx),
                                                                  ": ",
                                                                  "Added species data - '",
                                                                  paste(new_cols_chr, collapse = "', '"),
                                                                  "' - from '",
                                                                  spp_info_chr,
                                                                  "'.")
  # %%%

  obj$checkself()
  
  
  return(obj)
}

#' Joins external site information to a database object
#' 
#' @description 
#' This function joins external information to a [`Database`] object. Such 
#' information could be for example additional site-specific data that are not 
#' already present in the database.
#' @param obj A [`Database`] object.
#' @param site_info A [`dataframe`] or similar of categorical or numeric information 
#' to be added to the database.
#' @keywords database, external
#' @export
hib_join_sites <- function(obj, site_info){
  # Assertions
  assertthat::assert_that(
    inherits(obj, "Database"),
    is.data.frame(site_info)
    )

  # Tests
  # columns in site_info 
  # - at least 1 must be in obj$sites
  # - at least 1 must not be in obj$sites
  obj_site_cols <- colnames(obj$sites)
  info_cols <- colnames(site_info)
  
  n_matching_cols <- sum(info_cols %in% obj_site_cols)
  
  assertthat::assert_that(
    n_matching_cols > 0 & n_matching_cols < length(info_cols),
    msg = "Columns in site_info are not compatible with the current site data."
  )
  
  message(crayon::green("Joining site-related data into hib database...\n"))
  
  obj$sites <- dplyr::left_join(obj$sites,
                         site_info)

  obj$last_update <- lubridate::today()
  

  
  # %%%
  # TODO: would be useful to flag addition of data - source name/an id..
  site_info_chr <- deparse(substitute(site_info))
  new_cols_chr <- info_cols[!info_cols %in% obj_site_cols]

  # Append filtering flag to database object
  join_idx <- length(obj$processing_flags$joined_information) + 1
  
  obj$processing_flags$joined_information[[join_idx]] <- paste0(as.character(join_idx),
                                                                  ": ",
                                                                  "Added site data - '",
                                                                  paste(new_cols_chr, collapse = "', '"),
                                                                  "' - from '",
                                                                  site_info_chr, 
                                                                  "'.")
  # %%%

  obj$checkself()

  return(obj)
}