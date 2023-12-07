#' Save a hib object to RDS
#'
#' @description
#' The function saves the hibiiscus object in the specified file location.
#' @details
#' Saves a hibiiscus data object to file.
#' @param obj A hibiiscus object one of: `HibDatabase`, `HibLULegend`, 
#' `HibModelLME4`, `HibModelINLA`, `HibModelBRMS`, `HibCV`, `HibProjection`.
#' @param path A [`character`] pointing to a file (file names should have 
#' a ".rds" extension. Non-rds extensions will be amended to rds).
#' @param overwrite A [`logical`] indicating whether a previously saved file should
#' be overwritten. (Default: \code{FALSE})
#' @keywords data, save
#' @examples
#' \dontrun{
#'   hib_save(model_obj, "an/example/file_path/hib_model.rds", T)
#' }
#' @export
hib_save <- function(obj, path, overwrite = FALSE){
  # Add some checks...
  assertthat::assert_that(
    # Possible data types/classes of obj
    (inherits(obj, "HibDatabase") || 
      inherits(obj, "HibLULegend") || 
      inherits(obj, "HibModelLME4") ||
      inherits(obj, "HibModelINLA") || 
      inherits(obj, "HibModelBRMS") || 
      inherits(obj, "HibCV") ||
      inherits(obj, "HibProjection")) ,
    is.character(path),
    is.logical(overwrite)
  )
  
  obj$checkself()
  
  dir_path <- dirname(path)
  if (!dir.exists(dir_path)){
    message(crayon::yellow(paste0("\nWarning: '", dir_path, "' doesn't exist.",
                                  "\n",
                                  "This folder will be created.\n")))
    dir.create(dir_path, recursive = T)
  }
  
  if(!assertthat::has_extension(path, "rds")){
    message(crayon::yellow(paste0("\nWarning: specified file does not have an 'rds' extension.",
                                  "\n",
                                  "An 'rds' extension will be used for saving.\n")))
    ext <- tools::file_ext(path)
    if (ext != ""){
      path <- gsub(ext, "rds", path)
    } else {
      path <- paste0(path, ".rds")
    }
  }
  
  if (file.exists(path)){
    if (overwrite == FALSE){
      message(crayon::red(paste0("\nError: '", path, "' already exists.\n
                                    Data will NOT be saved.\n
                                    Use 'overwrite = TRUE' if you want to overwrite this file.\n")))
      stop("Hibiiscus data could not be saved. See above for details.", call. = F)
    }
    
    if (overwrite == TRUE){
      message(crayon::yellow(paste0("\nWarning: '", path, "' already exists.\n
                                    This file will be updated.\n")))
    }
  }
  
  
  message(crayon::green(paste0("\nSaving ", class(obj), " to: '", path, "' ...")) )
  saveRDS(obj, path)
  assertthat::assert_that(file.exists(path),
                          msg = "Writing a copy of the database failed...")
  
  invisible()
}





#' Load a hib object from RDS
#'
#' @description
#' The function loads the hibiiscus object fomr the specified file location.
#' @details
#' Loads a hibiiscus data object to file.
#' @param path A [`character`] pointing to a file (file names must have 
#' a ".rds" extension).
#' @keywords data, load
#' @examples
#' \dontrun{
#'   hib_obj <- hib_load("an/example/file_path/hib_model.rds")
#' }
#' @export
hib_load <- function(path){

  assertthat::assert_that(
    is.character(path),
    file.exists(path),
    assertthat::has_extension(path, "rds")
  )
  # Try and read in the database
  message(crayon::green(paste0("Loading Hibiiscus object from: '", path, "'...")))
  
  obj <- readr::read_rds(path)
  
  # Make some sensitivity checks that the object is correct
  assertthat::assert_that(
    # Possible data types/classes of obj
    (inherits(obj, "HibDatabase") || 
      inherits(obj, "HibLULegend") || 
      inherits(obj, "HibModelLME4") ||
      inherits(obj, "HibModelINLA") || 
      inherits(obj, "HibModelBRMS") || 
      inherits(obj, "HibCV") ||
      inherits(obj, "HibProjection"))
  )

  message(crayon::green(paste0(class(obj), " loaded...\n")) )
  
  obj$checkself()
  
  return(obj)
  
}