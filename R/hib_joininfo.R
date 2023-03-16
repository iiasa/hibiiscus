#' Joins external species information to a database object
#' 
#' @description 
#' This function joins external information to a [`Database`] object. Such information could be 
#' for example prepared data on species traits.
#' @param obj A [`Database`] object.
#' @param spp_info A [`dataframe`] or similar of categorical or numeric information to be added to the database.
#' @keywords database, external
hib_join_species <- function(obj, spp_info){
  # stop("TBD")
  obj$species <- left_join(obj$species,
                           spp_info)
  
  return(obj)
}

#' Joins external site information to a database object
#' 
#' @description 
#' This function joins external information to a [`Database`] object. Such information could be 
#' for example additional site-specific data that are not already present in the database.
#' @param obj A [`Database`] object.
#' @param site_info A [`dataframe`] or similar of categorical or numeric information to be added to the database.
#' @keywords database, external
hib_join_sites <- function(obj, site_info){
  # stop("TBD")
  obj$sites <- left_join(obj$sites,
                         site_info)
  
  return(obj)
}