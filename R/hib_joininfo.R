#' Joins external species information to a database object
#' 
#' @description 
#' This function joins external information to a [`Database`] object. Such information could be 
#' for example prepared data on species traits.
#' @param obj A [`Database`] object.
#' @param ids [`character`] giving the id column to be added.
#' @param info A [`vector`] of categorical or numeric information to be added to the database.
#' @keywords database, external
hib_join_species <- function(obj, ids, info){
  stop("TBD")
}

#' Joins external site information to a database object
#' 
#' @description 
#' This function joins external information to a [`Database`] object. Such information could be 
#' for example additional filtering data for sites not already present in the database.
#' @param obj A [`Database`] object.
#' @param ids [`character`] giving the id column to be added.
#' @param info A [`vector`] of categorical or numeric information to be added to the database.
#' @keywords database, external
hib_join_sites <- function(obj, ids, info){
  stop("TBD")
}