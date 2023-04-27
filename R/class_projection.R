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
                           is.data.frame(proj) || inherits(proj, "SpatRaster") || (is.list(proj) & !is.data.frame(proj)), # or a list of SpatRasters
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
                         # stop("Plotting functionality in development")
                         
                         # Get projection info - self$projection
                         proj <- self$projection
                         
                         # Check for self$projection_specs$Rescaled_prediction
                         resp_var <- "prediction"
                         resp_lab <- "Prediction"
                         if (self$projection_specs$Rescaled_prediction){
                           resp_var <- "rescaled_prediction"
                           resp_lab <- "Rescaled\nprediction"
                           
                         }
                         # If not rescaled..., just use a basic colour gradient
                         
                         # Initially, basic ggplot
                        
                         # TODO: add other plotting options...
                         
                         # TODO: inla projections have _mean etc...
                         print(inherits(self$model, "hib_inla_model"))
                         if (inherits(self$model, "hib_inla_model")){
                           resp_var <- paste(resp_var, "mean", sep = "_")
                         }
                         if (inherits(proj, "SpatRaster") ){
                           # terra::plot(proj$`2010`[[resp_var]])
                           proj <- terra::as.data.frame(proj[[resp_var]], xy = T)
                         }
                         if ((is.list(proj) & !is.data.frame(proj))){
                           proj <- purrr::map(proj,
                                              function(z){terra::as.data.frame(z[[resp_var]], xy = T)}) %>%
                             bind_rows(.id = "Year")
                         }
                         # TODO: make use of tidyterra??

                         plt <- ggplot(proj) +
                           geom_tile(aes(x = x, y = y,
                                         fill = !!sym(resp_var),
                                         colour = !!sym(resp_var))) +
                           labs(fill = resp_lab,
                                colour = resp_lab) +
                           theme_bw() +
                           theme(axis.text = element_text(size = 16),
                                 axis.title = element_text(size = 20),
                                 plot.title = element_text(size = 20),
                                 plot.subtitle = element_text(size = 18),
                                 strip.text.x = element_text(size = 18),
                                 strip.text.y = element_text(size = 18),
                                 strip.background = element_blank(),
                                 legend.text = element_text(size = 16),
                                 legend.title = element_text(size = 18),
                                 legend.text.align = 0)
                         
                         # If rescaled, can also improve colour scheme 
                         # 1 as midpoint, lower "bad" (0), higher "better" (max)
                         if (self$projection_specs$Rescaled_prediction){
                           max_resp <- max(proj[,resp_var])
                           
                           plt <- plt +
                             scale_fill_gradient2(low = "firebrick3",
                                                 mid = "gainsboro",
                                                 high = "dodgerblue3",
                                                 midpoint = 1,
                                                 limits = c(0,max_resp)) +
                             scale_colour_gradient2(low = "firebrick3",
                                                    mid = "gainsboro",
                                                    high = "dodgerblue3",
                                                    midpoint = 1,
                                                    limits = c(0,max_resp))
                             
                         } else if (self$projection_specs$Rescaled_prediction == F){
                           plt <- plt +
                             scale_fill_gradient(low = "gainsboro",
                                                  high = "dodgerblue3") +
                             scale_colour_gradient(low = "gainsboro",
                                                    high = "dodgerblue3")
                         }
                         
                         if ("Year" %in% colnames(proj)){
                           plt <- plt +
                             facet_wrap(~Year)
                         }
                         
                         plot(plt)
                         
                       }
                       
                     ),
                     # Private entries to be created
                     private = list(
                       siteids = list(),
                       length = list()
                     )
)


