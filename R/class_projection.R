# NOTE: really, just want to add projection data to model... 
# Sensible in terms of reproducibility...
# but seems inefficient in terms of memory...

#' Template object for a projection object
#' 
#' @description 
#' This is the template object for a projection object - this is the result of 
#' combining a fitted hib_model with data of interest to make a projection of 
#' biodiversity responses
#' @keywords internal
#' @import R6
#' @noRd
hib_mod_projection <- R6Class("hib_mod_projection",
                     # Public data objects and functions
                     public = list(
                       # version = NULL,
                       last_update = NULL,
                       model = NULL, # should this be a whole model class??
                       # input_data = NULL,
                       projection = NULL ,
                       projection_specs = NULL,
                       # Other things to add
                       # aoa, metric of dissimlarity
                       
                       
                       # processing_flags = NULL,
                       
                       # Executed when created
                       initialize = function(obj,
                                             proj,
                                             backtransform,
                                             rescale_pred){
                         # Check that input is either a character filepath
                         # Or dataframe/tibble
                         assertthat::assert_that(
                           inherits(obj, "hib_lme4_model") || inherits(obj, "hib_inla_model"), # obj is a hib model
                           is.data.frame(proj) || inherits(proj, "SpatRaster"), # or a list of SpatRasters
                           is.function(backtransform) || is.null(backtransform),
                           is.logical(rescale_pred)
                         )
                         
                         # ------ #
                         message(crayon::green("Initializing model projection class") )
                         
                         # Store projection object (list of data)
                         self$last_update <- lubridate::today()
                         
                         self$model <- obj
                         self$projection <- proj
                         self$projection_specs <- list("Backransformation_function" = backtransform,
                                                       "Rescaled_prediction" = rescale_pred)
                         
                         
                         # And check 
                         # Outsourced security checks
                         # self$checkself()
                         
                         # Add some data processing flags
                         # self$processing_flags <- list(measures_effort_corrected = F,
                         #                               abundance_rescaled = F,
                         #                               filtering_steps = list())
                       },
                       # Self check function that everything is in order
                       checkself = function(){
                         # ------ #
                         # Security checks - recheck
                         # assertthat::assert_that(
                         # )
                         # assertthat::assert_that(
                         # )
                         # invisible()
                       },
                       # Base printout function
                       print = function() {
                         stop("Print functionality in progress...")
                         # Compute some stats
                         # stats <- list(number_sites = length( unique(self$lu_df$SSBS) ),
                         #               number_lu = length( self$land_use ),
                         #               number_i = length( self$management ),
                         #               number_lui = length( self$land_use_management )
                         # )  
                         # 
                         # message(
                         #   paste0("Hibiiscus projection \n",
                         #          "(Last updated: ", crayon::green( self$last_update ), ")",
                         #          "\n",
                         #          "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                         #          "\t N land use types:              ", crayon::green(stats$number_lu), "\n",
                         #          "\t N management types:            ", crayon::green(stats$number_i), "\n",
                         #          "\t N land use x management types: ", crayon::green(stats$number_lui), "\n",
                         #          "\n"
                         #   )
                         # )
                         
                       },
                       # Save object
                       save_projection = function(output_path = NULL){
                         assertthat::assert_that(
                           is.character(output_path)
                         )
                         
                         if(!assertthat::has_extension(output_path, "rds")){
                           output_path <- paste0(output_path, ".rds")
                         }
                         
                         if (file.exists(output_path)){
                           message(crayon::yellow(paste0("Warning: '", output_path, "' already exists.",
                                                         "\n",
                                                         "This file will be updated.\n")))
                         }
                         
                         message(crayon::green(paste0("Saving projection to: '", output_path, "' ...\n")) )
                         
                         saveRDS(self, output_path)
                         assertthat::assert_that(file.exists(output_path),
                                                 msg = "Writing a copy of the projection failed...")
                         invisible()
                       },
                       
                       # TODO:
                       # summary plots - distr of luis... - either lu, i, or lui
                       plot = function(){
                         stop("Plotting functionality in development")
                         
                         # Get projection info - self$projection
                         # Check for self$projection$predicted_rescaled
                         
                         # Initially, basic ggplot
                       }
                       
                     ),
                     # Private entries to be created
                     private = list(
                       siteids = list(),
                       length = list()
                     )
)


