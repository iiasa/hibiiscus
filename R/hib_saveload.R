#' hib_save: Save Hib___ objects.
#'
#' @description
#' Saves a hibiiscus data object to disk.
#' 
#' @details
#' The function saves the hibiiscus object in the specified file location.
#' If the object (`HibProjection`) includes \code{SpatRaster} data, this will 
#' be separately saved in '.tif' file/s, with the same base-path.
#' 
#' @param obj A hibiiscus object, one of: `HibDatabase`, `HibLULegend`, 
#' `HibModelLME4`, `HibModelINLA`, `HibModelBRMS`, `HibCV`, `HibProjection`.
#' @param path A `character` pointing to a file (file names should have 
#' a ".rds" extension. Non-rds extensions will be amended to rds).
#' @param overwrite A `logical` indicating whether a previously saved file should
#' be overwritten. (Default: \code{FALSE})
#' 
#' @return NULL
#' 
#' @keywords data, models, save
#' 
#' @seealso
#' [hib_load()]
#' 
#' @examples
#' \dontrun{
#'   hib_save(model_obj, "an/example/file_path/hib_model.rds", T)
#' }
#' @export
hib_save <- function(obj, path, overwrite = FALSE){
  
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
  
  if ("checkself" %in% names(obj)){
    obj$checkself()
  }
  
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
  
  # print(file.exists(path))
  
  if (file.exists(path)){
    if (overwrite == FALSE){
      message(crayon::red(paste0("\nError: '", path, "' already exists.\n",
                                    "\tData will NOT be saved.\n",
                                    "Use 'overwrite = TRUE' if you want to overwrite this file.\n")))
      stop("Hibiiscus data could not be saved. See above for details.", call. = F)
    }
    
    if (overwrite == TRUE){
      message(crayon::yellow(paste0("\nWarning: '", path, "' already exists.\n",
                                    "\tThis file will be updated.\n")))
    }
  }
  
  # Note: Hib objects have two class items : the HibXXX, and R6, only want the first
  obj_class <- class(obj)
  obj_class <- obj_class[obj_class != "R6"]

  dc_obj <- NULL
  
  # check for rast (could be a list of rasters!) (only in HibProjection, $projection and $aoa)
  if (inherits(obj, "HibProjection")){
    # message("Checking HibProjection object for SpatRast data...")
    # If present, save to tif/s!
    
    # %%%
    # $projection 
    if (inherits(obj$projection, "SpatRaster")){
      
      # Deep-copy obj
      dc_obj <- obj$clone(deep = T)
      
      # specify output fp for raster
      proj_fp <- tools::file_path_sans_ext(path)
      proj_fp <- paste0(proj_fp, "_proj.tif")
      
      # store rast fp/s in new obj
      # dc_obj$projection_specs$projection_rast_fp <- proj_fp
      # TODO: incorp data-type!
      # save raster
      message("Saving HibProjection SpatRast data...")
      terra::writeRaster(obj$projection,
                         proj_fp,
                         overwrite = overwrite)
      # Dynamically store fps (basename only!)
      dc_obj$projection <- basename(proj_fp)
      
      # Save AOA here too - same format as projection
      if (!is.null(obj$aoa)) {
        aoa_fp <- tools::file_path_sans_ext(path)
        aoa_fp <- paste0(aoa_fp, "_aoa.tif")
        
        # store rast fp/s in new obj
        # dc_obj$aoa$aoa_rast_fp <- aoa_fp
        # TODO: incorp data-type!
        # save raster
        message("Saving HibProjection AOA SpatRast data...")
        terra::writeRaster(obj$aoa$proj_aoa,
                           aoa_fp,
                           overwrite = overwrite)
        # Dynamically store fps
        dc_obj$aoa$aoa$proj_aoa <- basename(aoa_fp)
      }
    }
    
    if (
      ( is.list(obj$projection) & 
        !is.data.frame(obj$projection) & 
        all(unlist(lapply(obj$projection, function(x){inherits(x, "SpatRaster")}))) )
    ){
      
      # Deep-copy obj
      dc_obj <- obj$clone(deep = T)
      
      # specify output fps for raster
      proj_fp <- tools::file_path_sans_ext(path)
      
      # Add in idx based on length of projection
      suff <- paste0("_proj_",
                     1:length(obj$projection),
                     ".tif")
      
      proj_fp <- paste0(proj_fp, suff)
      
      # store rast fp/s in new obj
      # dc_obj$projection_specs$projection_rast_fp <- proj_fp
      
      message("Saving HibProjection Projection SpatRaster data...")
      for (i in 1:length(obj$projection)){
        # TODO: incorp data-type!
        # save raster
        terra::writeRaster(obj$projection[[i]],
                           proj_fp[i],
                           overwrite = overwrite)
        # Dynamically store fps
        dc_obj$projection[[i]] <- basename(proj_fp[i])
      }
      
      # Save AOA here too - same format as projection
      if (!is.null(obj$aoa)) {
        # specify output fps for raster
        aoa_fp <- tools::file_path_sans_ext(path)
        
        # Add in idx based on length of projection
        suff <- paste0("_aoa_",
                       1:length(obj$aoa$proj_aoa),
                       ".tif")
        
        aoa_fp <- paste0(aoa_fp, suff)
        
        # store rast fp/s in new obj
        # dc_obj$aoa$aoa_rast_fp <- aoa_fp
        
        message("Saving HibProjection AOA SpatRaster data...")
        for (i in 1:length(obj$aoa$proj_aoa)){
          # TODO: incorp data-type!
          # save raster
          terra::writeRaster(obj$aoa$proj_aoa[[i]],
                             aoa_fp[i],
                             overwrite = overwrite)
          # Dynamically store fps
          dc_obj$aoa$proj_aoa[[i]] <- basename(aoa_fp[i])
        }
      }
    }
    
    # %%%
  }
  
  message(crayon::green(paste0("\nSaving ", obj_class, " to: '", path, "' ...")) )
  
  if (!is.null(dc_obj)) {
    
    saveRDS(dc_obj, path)
    
  } else {
    saveRDS(obj, path)  
  }
  

  assertthat::assert_that(file.exists(path),
                          msg = "Writing a copy of the database failed...")
  
  invisible()
}





#' hib_load: Load Hib___ objects.
#'
#' @description
#' Load a hibiiscus object from a saved RDS.
#' 
#' @details
#' The function loads the hibiiscus object from the specified file location.
#' If the original object (HibProjection) included \code{SpatRaster} data, this 
#' will be loaded from the appropriate tif file/s.
#' 
#' @param path A `character` pointing to a file (file names should have 
#' a ".rds" extension, if they don't, the extension will be added).
#' 
#' @return A hibiiscus object, one of: `HibDatabase`, `HibLULegend`, 
#' `HibModelLME4`, `HibModelINLA`, `HibModelBRMS`, `HibCV`, `HibProjection`
#' 
#' @keywords data, load
#' 
#' @seealso
#' [hib_save()]
#' 
#' @examples
#' \dontrun{
#'   hib_obj <- hib_load("an/example/file_path/hib_model.rds")
#' }
#' @export
hib_load <- function(path){

  assertthat::assert_that(
    is.character(path),
    file.exists(path)
    # assertthat::has_extension(path, "rds")
  )
  
  if (!assertthat::has_extension(path, "rds")){
    message(crayon::yellow(paste0("\nWarning: specified file does not have an 'rds' extension.",
                                  "\n",
                                  "An 'rds' extension will be added.\n")))
    ext <- tools::file_ext(path)
    if (ext != ""){
      path <- gsub(ext, "rds", path)
    } else {
      path <- paste0(path, ".rds")
    }
  }
  
  if (!file.exists(path)) {
    message(crayon::red(paste0("\nError: '", path, "' does not exist.\n",
                               "\tData will NOT be loaded.\n",)))
    stop("Hibiiscus data could not be loaded. See above for details.", call. = F)
  }
  
  # Try and read in the database
  message(crayon::green(paste0("\nLoading Hibiiscus object from: '", path, "'...\n")))
  
  obj <- readRDS(path)
  
  # Make some sensitivity checks
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

  obj_class <- class(obj)
  obj_class <- obj_class[obj_class != "R6"]
  
  message(crayon::green(paste0(obj_class, " loaded...\n")) )
  
  if ("checkself" %in% names(obj)){
    obj$checkself()
  }
  
  
  # check for rast data (could be a list of rasters!), only in HibProjection, ($projection, $aoa)
  # if present, aoa 'projection; should be self-similar to actual projection - dimensions, n years etc
  
  if (inherits(obj, "HibProjection")) {
    if (is.character(obj$projection)) {
      # print("here")
      message("Adding Projection SpatRaster ...")
      # Load proj basename
      # Add to path (minus path basename)
      rast_path <- file.path(dirname(path),
                             obj$projection)
      assertthat::assert_that(
        file.exists(rast_path)
      )
      
      # Put loaded rast in $projection slot
      obj$projection <- terra::rast(rast_path)
      
      assertthat::assert_that(
        inherits(obj$projection, "SpatRaster")
      )
      
      
      # If !is.null(aoa)
      if (!is.null(obj$aoa)) {
        message("Adding AOA SpatRaster ...")
        # Load aoa basename
        # Add to path (minus path basename)
        rast_path <- file.path(dirname(path),
                               obj$aoa$proj_aoa)
        
        assertthat::assert_that(
          file.exists(rast_path)
        )
        
        # Put loaded rast in $aoa$proj_aoa slot
        obj$aoa$proj_aoa <- terra::rast(rast_path)
        
        assertthat::assert_that(
          inherits(obj$aoa$proj_aoa, "SpatRaster")
        )
        
      }
      
    }
    
    if (
      ( is.list(obj$projection) & 
        !is.data.frame(obj$projection) & 
        all(unlist(lapply(obj$projection, function(x){inherits(x, "character")}))) )
    ) {
      print("here, list")
      
      message("Adding Projection SpatRasters ...")
      
      for (i in 1:length(obj$projection)) {
        # Load proj basename
        # Add to path (minus path basename)
        rast_path <- file.path(dirname(path),
                               obj$projection[[i]])
        assertthat::assert_that(
          file.exists(rast_path)
        )
        
        # Put loaded rast in $projection slot
        obj$projection[[i]] <- terra::rast(rast_path)
        
        assertthat::assert_that(
          inherits(obj$projection[[i]], "SpatRaster")
        )
      }
      
      
      # If !is.null(aoa)
      if (!is.null(obj$aoa)) {
        message("Adding AOA SpatRasters ...")
        
        for (i in 1:length(obj$aoa$proj_aoa)) {
          # Load aoa basename
          # Add to path (minus path basename)
          rast_path <- file.path(dirname(path),
                                 obj$aoa$proj_aoa[[i]])
          
          assertthat::assert_that(
            file.exists(rast_path)
          )
          
          # Put loaded rast in $aoa$proj_aoa slot
          obj$aoa$proj_aoa[[i]] <- terra::rast(rast_path)
          
          assertthat::assert_that(
            inherits(obj$aoa$proj_aoa[[i]], "SpatRaster")
          )
        }
      }
      
    }
  }
  
  return(obj)
  
}