#' Create an lu legend
#' 
#' @description 
#' This function creates a [`R6`] object containing the loaded land-use legend. 
#' 
#' This object is then used in other functions of this package to prepare or filter
#' associated databases and enable ecological impact analyses. 
#' @param input Either a [`character`] pointing to a saved RDS file, or a [`dataframe`]-like
#' object within the environment.
#' @keywords database, land-use, legend
#' @examples
#' \dontrun{
#'   lu_leg <- create_lu_legend("../Data/landuse_dataframe.RDS")
#'   lu_leg
#' }
#' @export
create_lu_legend <- function(input){
  # %%%
  assertthat::assert_that(
    is.character(input) || is.data.frame(input)
  )
  if (is.character(input)){
    assertthat::assert_that(file.exists(input))
  }
  
  # ------ #
  message("Creating land-use legend...\n")
  
  leg <- lu_legend$new(input = input)
  
  return(leg)
  # %%%
  
}

#' Save a previously created legend
#' 
#' @description 
#' The function [`save_legend`] saves the output of the formatted [`legend`] object in a specific
#' place.
#' @details 
#' Expecting a file with \code{rds} ending.
#' @param leg A land-use [`legend`].
#' @param output_path A [`character`] pointing to a filename.
#' @keywords save, land-use, legend
#' @export
save_legend <- function(leg, output_path) leg$save_legend(output_path)



#' Load a saved lu legend
#' 
#' @description 
#' This function loads a previously created [`R6`] object containing a land-use legend class. 
#' 
#' @param input A [`character`] pointing to a saved RDS file.
#' @keywords database, land-use, legend
#' @examples
#' \dontrun{
#'   lu_leg <- load_lu_legend("../Data/landuse_legend.RDS")
#'   lu_leg
#' }
#' @export
load_lu_legend <- function(path_input){
  assertthat::assert_that(
    is.character(path_input),
    file.exists(path_input),
    assertthat::has_extension(path_input, "rds")
  )
  # Try and read in the database
  message(crayon::green(paste0("Loading land-use legend from: '", path_input, "...\n")))
  
  leg <- readr::read_rds(path_input)
  
  # TODO: Make some sensitivity checks that the object is correct
  leg$checkself()
  
  return(leg)
  
}


# TODO:
# Need to incl function to add legend...
# Add to a list in database?


#' Set a loaded lu legend as the LU legend within a hibiiscus database
#'
#' @description
#' This function adds a land-use legend to a database - specifically modifying the 
#' sites dataframe.
#'
#' @param obj A [`hib_database`].
#' @param leg A [`hib_lu_legend`].
#' @keywords database, land-use, legend
#' @examples
#' \dontrun{
#'   obj <- set_lu_legend(obj, leg)
#'   obj
#' }
#' @export
set_lu_legend <- function(obj, leg){
  assertthat::assert_that(
    inherits(obj, "Database") ,
    inherits(leg, "lu_legend") 
  )
  # Try and read in the database
  message(crayon::green("Setting land-use legend for hib database...\n"))
  
  # Make some sensitivity checks that the inputs are okay
  obj$checkself()
  leg$checkself()
  
  # check for column names already in sites (Predominant_land_use, Use_intensity, LUI)
  cols_to_drop <- colnames(obj$sites)[colnames(obj$sites) %in% c("Predominant_land_use", "Use_intensity", "LUI")]
  # Drop cols from sites
  obj$sites <- obj$sites %>%
    select(!all_of(cols_to_drop)) %>%
    # left join new legend info
    dplyr::left_join(leg$lu_df)
  
  return(obj)
  
}


