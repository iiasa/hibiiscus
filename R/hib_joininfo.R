#' hib_join_species: Join external species information to a HibDatabase object
#' 
#' @description 
#' This function joins external information to a [`HibDatabase`] object. Such 
#' information could be for example prepared data on species traits.
#' @param obj A `HibDatabase` object.
#' @param spp_info A `data.frame` or similar of categorical or numeric information 
#' to be added to the species data in the `HibDatabase`. 
#' Must share at least one column with `HibDatabase$species`, e.g. 'Best_guess_binomial'.
#' 
#' @return HibDatabase
#' 
#' @keywords database, external
#' 
#' @examples 
#' \dontrun{
#'   hib_db <- hib_join_species(hib_db, new_spp_data_df)
#' } 
#' @export
hib_join_species <- function(obj, spp_info){

  # Assertions
  assertthat::assert_that(
    inherits(obj, "HibDatabase"),
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
  message(crayon::green("Joining species-related data into HibDatabase...\n"))
  
  obj$species <- dplyr::left_join(obj$species,
                                  spp_info)  
  
  obj$last_update <- lubridate::now()
  
  # %%%
  # Flag addition of data - source name/an id..
  spp_info_chr <- deparse(substitute(spp_info))
  new_cols_chr <- info_cols[!info_cols %in% obj_spp_cols]

  # Append join flag to database object
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


#' hib_join_sites: Join external site information to a HibDatabase object
#' 
#' @description 
#' This function joins external information to a [`HibDatabase`] object. Such 
#' information could be for example additional site-specific data that are not 
#' already present in the database.
#' 
#' @param obj A `HibDatabase` object.
#' @param site_info A `data.frame` or similar of categorical or numeric information 
#' to be added to the site dataframe within the `HibDatabase`. 
#' Must share at least one column with `HibDatabase$sites`, e.g., 'SSBS'.
#' 
#' @return HibDatabase
#' 
#' @keywords database, external
#' 
#' @examples 
#' \dontrun{
#'   hib_db <- hib_join_sites(hib_db, new_site_data_df)
#' } 
#' @export
hib_join_sites <- function(obj, site_info){
  # Assertions
  assertthat::assert_that(
    inherits(obj, "HibDatabase"),
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
  
  message(crayon::green("Joining site-related data into HibDatabase...\n"))
  
  obj$sites <- dplyr::left_join(obj$sites,
                         site_info)

  obj$last_update <- lubridate::now()
  
  # %%%
  # Flag addition of data - source name/an id..
  site_info_chr <- deparse(substitute(site_info))
  new_cols_chr <- info_cols[!info_cols %in% obj_site_cols]

  # Append join flag to database object
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


#' hib_join_site_contrasts: Join external site-contrast information to a HibDatabase object
#' 
#' @description 
#' This function joins external information to a [`HibDatabase`] object. Such 
#' information could be for example additional site-contrast-specific data that 
#' are not already present in the database.
#' 
#' @param obj A `HibDatabase` object.
#' @param site_contrast_info A `data.frame` or similar of categorical or numeric 
#' information to be added to the database. 
#' Must share at least one column with `HibDatabase$site_contrasts`. 
#' Ideally, both/one of 'SSBS1' and 'SSBS2' are present.
#' 
#' @return HibDatabase
#' 
#' @keywords database, external
#' 
#' @examples 
#' \dontrun{
#'   hib_db <- hib_join_site_contrasts(hib_db, new_site_contrast_data_df)
#' } 
#' @export
hib_join_site_contrasts <- function(obj, site_contrast_info){
  # Assertions
  assertthat::assert_that(
    inherits(obj, "HibDatabase"),
    is.data.frame(site_contrast_info)
    )

  # if site_contrasts exists, left_join
  if (!is.null(obj$site_contrasts)){
    # Tests
    # columns in site_contrast_info 
    # - at least 1 must be in obj$site_contrasts
    # - at least 1 must not be in obj$site_contrasts
    
    obj_site_contr_cols <- colnames(obj$site_contrasts)
    info_cols <- colnames(site_contrast_info)
    
    n_matching_cols <- sum(info_cols %in% obj_site_contr_cols)
    
    assertthat::assert_that(
      n_matching_cols > 0 & n_matching_cols < length(info_cols),
      msg = "Columns in site_contrast_info are not compatible with the current site_contrast data."
    )
    
    message(crayon::green("Joining site-contrast-related data into HibDatabase...\n"))
    
    obj$site_contrasts <- dplyr::left_join(obj$site_contrasts,
                                           site_contrast_info)
  }

  # If it doesn't, stop
  if (is.null(obj$site_contrasts)){
    # Check SSBS1, SSBS2...
    # Refuse to create
    message(crayon::red("Site-contrast data not present in HibDatabase, cannot join new data...\n"))
    stop("hib_join_site_contrasts unsuccessful. See above for details.", call. = F)
  }

  obj$last_update <- lubridate::now()
  
  # %%%
  # Flag addition of data - source name/an id..
  site_contr_info_chr <- deparse(substitute(site_contrast_info))
  new_cols_chr <- info_cols[!info_cols %in% obj_site_contr_cols]

  # Append join flag to database object
  join_idx <- length(obj$processing_flags$joined_information) + 1
  
  obj$processing_flags$joined_information[[join_idx]] <- paste0(as.character(join_idx),
                                                                  ": ",
                                                                  "Added site-contrast data - '",
                                                                  paste(new_cols_chr, collapse = "', '"),
                                                                  "' - from '",
                                                                  site_contr_info_chr, 
                                                                  "'.")
  # %%%

  obj$checkself()

  return(obj)
}
