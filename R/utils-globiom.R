#' prep_globiom: Prepare GLOBIOM output for projection.
#'
#' @description Function to format a prepared GLOBIOM netCDF file for use in 
#' Hibiiscus projections.
#' 
#' @details This function expects a downscaled GLOBIOM output. Variables are
#' extracted and converted to raster format, with each raster stack representing
#' a particular year.
#' Likely of little use for anyone outside IIASA.
#'
#' @param file_path A filepath (`character`) pointing to a GLOBIOM output in netCDF format.
#' @param out_as_rast A `logical` denoting if the output type should be a (list)
#' of `SpatRaster`, or a `stars` object (default: \code{TRUE}, i.e. returns  `SpatRaster`).
#' @param incl_pixel_area A [`logical`] denoting if the output should include the
#' 'pixel_area_ha' variable (Default: \code{TRUE}).
#' @param vars_to_include A `character` vector of variables to be extracted. 
#' If \code{NULL}, all variables will be extracted.
#' @param vars_to_ignore A `character` vector of variables to be ignored 
#' If \code{NULL}, all variables will be extracted - unless vars_to_include is not \code{NULL}.
#' 
#' @return A (list of) `SpatRaster` stacks with the formatted GLOBIOM predictors.
#'
#' @keywords utils
#'
#' @examples
#' \dontrun{
#' # Expects a filename pointing to a globiom-netCDF file.
#' globiom_rast <- prep_globiom("path/to/globiom_data.nc")
#' }
#'
#' @export
prep_globiom <- function(file_path,
                         out_as_rast = TRUE, # Currently this is enforced
                         incl_pixel_area = TRUE,
                         vars_to_include = NULL,
                         vars_to_ignore = NULL){
  # Assertions
  assertthat::assert_that(
    is.character(file_path),
    file.exists(file_path),
    assertthat::has_extension(file_path, "nc"),
    is.logical(out_as_rast),
    is.logical(incl_pixel_area),
    is.null(vars_to_include) || is.character(vars_to_include),
    is.null(vars_to_ignore) || is.character(vars_to_ignore)
  )
  
  # Get vars and dims
  tmp_nc <- ncdf4::nc_open(file_path)
  dims_ <- names(tmp_nc$dim)
  vars_ <- names(tmp_nc$var)
  
  if (!is.null(vars_to_include)){
    vars_ <- vars_[vars_ %in% vars_to_include]
  }
  if (!is.null(vars_to_ignore)){
    vars_ <- vars_[!vars_ %in% vars_to_ignore]
  }
  
  if (incl_pixel_area){
    # Check that "pixel_area_ha" is in vars_
    # If not, add
    
    if (!any(grepl("pixel_area", vars_, ignore.case = TRUE))){
      all_vars <- names(tmp_nc$var)
      pixel_area_var <- all_vars[grepl("pixel_area", all_vars, ignore.case = TRUE)]
      if (length(pixel_area_var) == 1){
        vars_ <- c(vars_, pixel_area_var)
      } 
      if (length(pixel_area_var) == 0){
        message(crayon::yellow("\nWarning: pixel area variable not found."))
        incl_pixel_area <- F
      }
      if (length(pixel_area_var) > 1){
        # Drop pixel_area[_ha]
        message(crayon::yellow(paste0("\nWarning: ",
                                      length(pixel_area_var), 
                                      " possible pixel area variables found.",
                                      "Neither will be included.")))
        incl_pixel_area <- F
      }
    }
  }
  
  if (!incl_pixel_area){
    # Drop pixel_area[_ha]
    vars_ <- vars_[!grepl("pixel_area", vars_, ignore.case = TRUE)]
  }
  
  assertthat::assert_that(
    length(vars_) >= 1,
    msg = "Variable selection has removed all variables."
  )  
    
  if (!grepl("netcdf", stars::detect.driver(file_path), ignore.case = TRUE)){
    stop("File format not recognized!")
  }
  
  message(crayon::green(paste0("\nLoading GLOBIOM NetCDF: ", file_path)))
  

  # For each vars_, load ?
  attr_ls <- list()
  stars_res <- vector()
  stars_area <- NULL
  
  for (v in vars_){
    # v <- vars_[1]
    # tmp_attr <- ncdf4::ncatt_get(tmp_nc, varid = v, verbose = FALSE)
    attr_ls[[v]] <- ncdf4::ncatt_get(tmp_nc, varid = v, verbose = FALSE)
    tmp_star <- stars::read_ncdf(file_path,
                                 var = v,
                                 proxy = FALSE,
                                 make_time = TRUE, # Make time on 'time' band
                                 make_units = FALSE) # To avoid unnecessary errors due to unknown units
    
    # if v is pixel_area...
    # store this separately...
    if (grepl("pixel_area", v, ignore.case = TRUE)){
      if(tmp_nc$var[[v]]$units %in% c("km2","ha","m2")){
        stars_area <- tmp_star[drop = TRUE]
      }
      rm(tmp_star)
      next()
    }
    
    # Get dimensions other than x, y and time 
    # Use to split stars object into layers...
    # Commonly used column names
    check <- c("x", "X", "lon", "longitude", 
              "y", "Y", "lat", "latitude", 
              "time", "Time", "year", "Year")
    extra_dims <- names(stars::st_dimensions(tmp_star))[!names(stars::st_dimensions(tmp_star)) %in% check]
    
    if(length(extra_dims) > 0){
      for(dim_nm in extra_dims){
        # FIXME: Dirty hack to remove forest zoning
        # if(length( grep("zone", dim_nm , ignore.case = T) ) > 0) next()
        
        # And class units as description from over
        dim_units <- tmp_nc$dim[[dim_nm]]$units
        dim_units <-  dim_units |>
          strsplit("/") |>
          # Remove emptyspace and, ?numbers? special symbols
          # 0-9
          sapply(function(y)  gsub("[^A-Za-z///' ]", "" , y, ignore.case = TRUE) ) |>
          sapply(function(y)  gsub(" ", "" , y, ignore.case = TRUE) )
        # Convert to vector and make names
        dim_units <- paste0(
          v, "__",
          make.names(unlist(dim_units)) |> as.vector()
        )
        
        tmp_star <- tmp_star |> stars:::split.stars(dim_nm) |> stats::setNames(nm = dim_units)
        
        tmp_dim_spec <- stars::st_dimensions(tmp_star)
        # Need to store 
        
        # FIXME: Dirty hack to deal with the forest zone dimension
        # If there are more dimensions than 3, aggregate over them
        # if( length(stars::st_dimensions(tmp_star)) >3){
        #   # Aggregate spatial-temporally
        #   tmp_star <- stars::st_apply(ff, c("longitude", "latitude", "time"), sum, na.rm = TRUE)
        # }
      }
    }
    
    # Add 
    stars_res <- c(stars_res, tmp_star)
    rm(tmp_star)
  }
  invisible(gc())
  
    # If only interested in getting stars_area, return just this, as rast...
  if (length(stars_res) == 0 &
      !is.null(stars_area)){
    assertthat::assert_that(
      length(dim(stars_area)) == 2
    )
    # stars_area <- stars_area[drop = T]
    rast_area <- terra::rast(stars_area)
    return(rast_area)
  }
  
  # Process stars results
  assertthat::assert_that(length(stars_res)>0)
  stars_res <- stars::st_as_stars(stars_res)


  # Need to set dimension names...
  # Just use tmp_dim_spec, check that length of tmp_dim_spec is equal to length of 
  if (length(stars::st_dimensions(stars_res)) == length(tmp_dim_spec)){
    stars::st_dimensions(stars_res) <- tmp_dim_spec # Target dimensions
  } else {
    stop("Star dimensions do not match!")
  }
  # Martin's approach : 
  # dim_spec <- tmp_dim_spec[c(
  #   grep("x|lon",names(tmp_dim_spec), ignore.case = TRUE,value = TRUE),
  #   grep("y|lat",names(tmp_dim_spec), ignore.case = TRUE,value = TRUE),
  #   grep("year|time",names(tmp_dim_spec), ignore.case = TRUE,value = TRUE)
  # )] # Order assumed to be correct
  
  assertthat::assert_that(length(names(stars::st_dimensions(stars_res)))==3)
  # ??
  
  # For each year, rasterise 
  # Return as a list of raster stacks - incl cell area (optional) - in-case weighting by this...
  out <- stars_to_rast(stars_res)
  
  if (incl_pixel_area){
    # Convert pixel area stars to raster
    assertthat::assert_that(
      length(dim(stars_area)) == 2
    )
    # stars_area <- stars_area[drop = T]
    rast_area <- terra::rast(stars_area)
    # bind to each item in out...
    out <- lapply(out, function(x, r = rast_area){c(x, r)})
  }
          
  return(out)
  
  # TODO: what would a stars output look like? - LU layers + lat/lon/time dimentions?
}


#' stars_to_rast: Converts a stars object to list of rasters
#'
#' @description This is a small helper function to convert a [`stars`] object to
#' a [`SpatRaster`] object. If multiple time periods are in the [`stars`] object,
#' then a [`list`] of [`SpatRaster`] (one item per time) will be returned.
#'
#' @param obj A [`stars`] object with a \code{"time"} dimension plus x/long and y/lat.
#'
#' @returns A [`list`] containing [`SpatRaster`] objects.
#'
#' @keywords scenario, internal
#'
#' @noRd
stars_to_rast <- function(obj){
  assertthat::assert_that(
    inherits(obj, 'stars')
  )
  
  n_dim <- length(dim(obj))
  
  # If no time dim - i.e., single timestep
  if (n_dim < 3){
    # If only one time point, return single raster stack object..
    out <- terra::rast(obj)
    return(out)
  }
  
  # If time, separate stars object by time
  # time_dim_idx <- which(names(dim(obj)) == "time")
  assertthat::assert_that(
    "time" %in% names(dim(obj))
  )
  
  time_dim_vals <- stars::st_get_dimension_values(obj, "time", center = TRUE)
  # stars::st_get_dimension_values(stars_res, "time", center = TRUE)
  time_idx <- 1:length(time_dim_vals)
  
  out <- vector("list", length(time_idx))
  names(out) <- as.character(time_dim_vals)
  
  for(t_ in time_idx){
    # Slice to a specific time frame for each
    tmp_rast <- obj |> stars:::slice.stars("time", t_) |>
      terra::rast() # Or alternatively rast
    names(tmp_rast) <- names(obj)
    
    # Set times for rast
    terra::time(tmp_rast) <- paste0(time_dim_vals[t_], "-01-01") |> 
      as.POSIXct(tz = "UTC") |> 
      as.Date() |> 
      rep(times = terra::nlyr(tmp_rast))
    
    out[[t_]] <- tmp_rast
  }
  
  return(out)
}

