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
                       # aoa, metric of dissimilarity
                       
                       
                       # processing_flags = NULL,
                       
                       # Executed when created
                       initialize = function(obj,
                                             proj,
                                             backtransform,
                                             rescale_pred){
                         # Check that input is either a character filepath
                         # Or dataframe/tibble
                         assertthat::assert_that(
                           inherits(obj, "hib_lme4_model") || inherits(obj, "hib_inla_model") || inherits(obj, "hib_brms_model"), # obj is a hib model
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
                         # stop("Print functionality in progress...")
                         # Compute some stats/get some info
                         
                         mod_type <- ""
                         pred_vars <- ""
                         n_pred_vars <- NA
                         if (inherits(self$model, "hib_lme4_model")){
                           mod_type <- "lme4"
                           
                           all_covars <- all.vars(formula(delete.response(terms(formula(self$model$model)))))
                           rand_vars <- names(lme4::ranef(self$model$model))
                           fix_vars <- all_covars[!all_covars %in% rand_vars]
                         }
                         
                         if (inherits(self$model, "hib_inla_model")){
                           mod_type <- "INLA"
                           
                           fix_vars <- all.vars(formula(delete.response(terms(self$model$user_formula$response_fixed))))
                           fix_vars <- fix_vars[!fix_vars == "Intercept"]
                           fix_vars <- fix_vars[!fix_vars == "(Intercept)"]
                           
                           if (self$model$inla_incl_spatial_term){
                             fix_vars <- c(fix_vars, "spatial-random-effect")
                           }
                         }
                         
                         if (inherits(self$model, "hib_brms_model")){
                           mod_type <- "brms"
                           # print(self$model$formula)
                           all_covars <- all.vars(formula(delete.response(terms(formula(self$model$model$formula)))))
                           rand_vars <- names(brms::ranef(self$model$model))
                           fix_vars <- all_covars[!all_covars %in% rand_vars]
                         }
                         
                         
                         n_pred_vars <- length(fix_vars)
                         if (length(fix_vars) > 3){
                           pred_vars <- paste(paste(fix_vars[1:3], collapse = ", ", ", ...", sep = ""))
                           pred_vars <- paste0(pred_vars, " (", n_pred_vars, " in total)")
                         } else {
                           pred_vars <- paste(fix_vars, collapse = ", ")
                         }
                         
                     
                         # spatial extent
                         # temporal extent
                         time_rng <- c(NA, NA)
                         # n_yrs <- NA
                         if (is.data.frame(self$projection)){
                           coords <- data.frame(x = NA, y = NA)
                           if (("Longitude" %in% colnames(self$projection)) &
                               ("Latitude" %in% colnames(self$projection)) ){
                             coords <- self$projection[,c("Longitude", "Latitude")]
                           }
                           if (("Longitude_midpoint" %in% colnames(self$projection)) &
                               ("Latitude_midpoint" %in% colnames(self$projection)) ){
                             coords <- self$projection[,c("Longitude_midpoint", "Latitude_midpoint")]
                           }
                           if (("x" %in% colnames(self$projection)) &
                               ("y" %in% colnames(self$projection)) ){
                             coords <- self$projection[,c("x", "y")]
                           }
                           if (("X" %in% colnames(self$projection)) &
                               ("Y" %in% colnames(self$projection)) ){
                             coords <- self$projection[,c("X", "Y")]
                           }
                           
                           x_rng <- range(coords[,1])
                           y_rng <- range(coords[,2])
                           
                           
                           if ("Year" %in% colnames(self$projection)){
                             time_rng <- range(self$projection$Year)
                           }
                           
                         }
                         if (inherits(self$projection, "SpatRaster")){
                           xy_rng <- terra::ext(self$projection)
                           x_rng <- xy_rng[1:2]
                           y_rng <- xy_rng[3:4]
                         }
                         if ((is.list(self$projection) & !is.data.frame(self$projection))){

                           xy_rng_df <- purrr::map_dfr(.x = self$projection,
                                          .f = function(.x){
                                            xy_rng <- terra::ext(.x)
                                            data.frame(xmin = xy_rng[1],
                                                       xmax = xy_rng[2],
                                                       ymin = xy_rng[3],
                                                       ymax = xy_rng[4])
                                          })
                           x_rng <- c(min(xy_rng_df$xmin), max(xy_rng_df$xmax))
                           y_rng <- c(min(xy_rng_df$ymin), max(xy_rng_df$ymax))
                           
                           time_rng <- range(as.integer(names(self$projection)))
                         }
                         
                         # n yrs
                         time_rng <- unique(time_rng)
                         if (length(time_rng) > 1){
                           n_yrs <- diff(time_rng) 
                         } else {
                           n_yrs <- 1
                         }
                         if (all(is.na(time_rng))){
                           n_yrs <- NA
                         }

                         # prnt_info <- list(mod_type = mod_type,
                         #                    )
                         message(
                           paste0("Hibiiscus projection \n",
                                  "(Last updated: ", crayon::green( self$last_update ), ")\n",
                                  "\n",
                                  "Model type:        ", mod_type, "\n",
                                  "\n",
                                  "Spatial extent:\n",
                                  "\t Longitude:      ", paste(x_rng, collapse = ", "), "\n",
                                  "\t Latitude:       ", paste(y_rng, collapse = ", "), "\n",
                                  "\n",
                                  "Temporal extent: \t ", paste(time_rng, collapse = "--"), " (", n_yrs, " year/s)", "\n",
                                  "\n",
                                  "Predictors:      \t ", pred_vars, "\n",
                                  "\n"
                           )
                         )
                         
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
                       
                       plot = function(type = "response-map"){
                         
                         type <- match.arg(type,
                                            c("response-map", "bivariate-map",
                                              "response-distr", "bivariate-distr"))

                         basic_thm <- ggplot2::theme(axis.text = ggplot2::element_text(size = 16),
                                           axis.title = ggplot2::element_text(size = 20),
                                           plot.title = ggplot2::element_text(size = 20),
                                           plot.subtitle = ggplot2::element_text(size = 18),
                                           strip.text.x = ggplot2::element_text(size = 18),
                                           strip.text.y = ggplot2::element_text(size = 18),
                                           strip.background = ggplot2::element_blank(),
                                           legend.text = ggplot2::element_text(size = 16),
                                           legend.title = ggplot2::element_text(size = 18),
                                           legend.text.align = 0) 
                         # TODO: add arg (type)
                         # Basic is response_map
                         # add bivariate_map - only currently works for INLA/BRMS
                         # add response_distr
                         # add bivariate_distr - only currently works for INLA/BRMS
                         
                         
                         # TODO: add bivariate plot - mean response v sd of response...
                         # in scaled space - 
                         
                         # INLA - (rescaled_)prediction_sd
                         # BRMS - (rescaled_)prediction_error
                         
                         # Get projection info - self$projection
                         proj <- self$projection
                         
                         # Check for self$projection_specs$Rescaled_prediction
                         resp_var <- "prediction"
                         resp_lab <- "Prediction"
                         if (self$projection_specs$Rescaled_prediction){
                           resp_var <- "rescaled_prediction"
                           resp_lab <- "Rescaled\nprediction"
                         }
                         if (inherits(self$model, "hib_inla_model")){
                           resp_var <- paste(resp_var, "mean", sep = "_")
                         }
                         
                         vars_ <- c(resp_var)
                         
                         if (type %in% c("bivariate-map", "bivariate-distr")){

                           if (inherits(self$model, "hib_lme4_model")){
                             message(crayon::red("\nBivariate plots not yet supported for projections from lme4-type models ...\n")) 
                             stop() 
                           }
                           if (inherits(self$model, "hib_inla_model")){
                             sd_var <- "prediction_sd"
                           }
                           if (inherits(self$model, "hib_brms_model")){
                             sd_var <- "prediction_error"
                           }
                           
                           if (self$projection_specs$Rescaled_prediction){
                             sd_var <- paste0("rescaled_", sd_var)
                           }
                           
                           if (type == "bivariate-map"){
                             resp_lab <- "Bivariate plot"  
                           }
                           if (type == "bivariate-distr"){
                             resp_lab <- gsub("\n", " ", resp_lab)
                             sd_lab <- "Uncertainty"
                           }
                           
                           
                           vars_ <- c(vars_, sd_var)
                           
                         }
                         
                         
                         # If not rescaled..., just use a basic colour gradient
                         # Initially, basic ggplot
                        
                         # TODO: add other plotting options...
                         # such as...
                         
                         if (is.data.frame(proj) ){
                           # terra::plot(proj$`2010`[[resp_var]])
                           vars__ <- c("x", "y", "Year", vars_)
                           vars__ <- vars__[vars__ %in% colnames(proj)]
                           proj <- proj[vars__]
                         }
                         if (inherits(proj, "SpatRaster") ){
                           # terra::plot(proj$`2010`[[resp_var]])
                           proj <- terra::as.data.frame(proj[[vars_]], xy = T)
                         }
                         if ((is.list(proj) & !is.data.frame(proj))){
                           proj <- purrr::map(proj,
                                              function(z){terra::as.data.frame(z[[vars_]], xy = T)}) %>%
                             dplyr::bind_rows(.id = "Year")
                         }
                         # TODO: make use of tidyterra??
                         # print(head(proj))

                         
                         if (type == "response-map"){
                           plt <- ggplot2::ggplot(proj) +
                             ggplot2::geom_tile(aes(x = x, y = y,
                                                    fill = !!sym(resp_var),
                                                    colour = !!sym(resp_var))) +
                             ggplot2::labs(fill = resp_lab,
                                           colour = resp_lab) +
                             ggplot2::coord_fixed() +
                             ggplot2::theme_bw() +
                             basic_thm
                           
                           # If rescaled, can also improve colour scheme 
                           # 1 as midpoint, lower "bad" (0), higher "better" (max)
                           if (self$projection_specs$Rescaled_prediction){
                             max_resp <- max(proj[,resp_var])
                             
                             plt <- plt +
                               ggplot2::scale_fill_gradient2(low = "firebrick3",
                                                             mid = "gainsboro",
                                                             high = "dodgerblue3",
                                                             midpoint = 1,
                                                             limits = c(0,max_resp)) +
                               ggplot2::scale_colour_gradient2(low = "firebrick3",
                                                               mid = "gainsboro",
                                                               high = "dodgerblue3",
                                                               midpoint = 1,
                                                               limits = c(0,max_resp))
                             
                           } else if (self$projection_specs$Rescaled_prediction == F){
                             plt <- plt +
                               ggplot2::scale_fill_gradient(low = "gainsboro",
                                                            high = "dodgerblue3") +
                               ggplot2::scale_colour_gradient(low = "gainsboro",
                                                              high = "dodgerblue3")
                           }
                           
                           if ("Year" %in% colnames(proj)){
                             plt <- plt +
                               ggplot2::facet_wrap(~Year)
                           }
                         }
                         
                         if (type == "bivariate-map"){
                           
                           # TODO: add colour scale options...
                           # c("Bluegill", "BlueGold", "BlueOr", "BlueYl", "Brown",
                           #   "Brown2", "DkBlue","DkBlue2", "DkCyan","DkCyan2", "DkViolet",
                           #   "DkViolet2", "GrPink","GrPink2", "PinkGrn", "PurpleGrn", "PurpleOr")
                           legend <- biscale::bi_legend(pal = "GrPink2",
                                                        dim = 3,
                                                        xlab = "Larger response",
                                                        ylab = "Larger uncert.",
                                                        size = 16)
                           
                           # Normalise response and sd...??
                           
                           proj$var1 <- proj[[resp_var]]
                           proj$var2 <- proj[[sd_var]]
                           
                           suppressWarnings(
                             proj <- biscale::bi_class(proj, x = var1, y = var2, dim = 3, style = "quantile")
                           )
                           
                           plt <- ggplot2::ggplot() +
                             ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = bi_class),
                                                data = proj) +
                             biscale::bi_theme(base_size = 16) +
                             biscale::bi_scale_fill(pal = "GrPink2", na.value = "transparent") +
                             # coord_quickmap() +
                             ggplot2::labs(
                               title = resp_lab,
                               x = "",
                               y = ""
                             ) +
                             ggplot2::theme_bw() +
                             basic_thm +
                             ggplot2::theme(legend.position = "none")
                           
                           if ("Year" %in% colnames(proj)){
                             plt <- plt +
                               ggplot2::facet_wrap(~Year)
                           }
                           
                           plt <- cowplot::ggdraw() +
                             cowplot::draw_plot(plt, 0, 0, 0.75, 1) +
                             cowplot::draw_plot(legend, 0.75, .45, 0.2, 0.2)
                         }
                         
                         # distribution of response... - histogram
                         if (type == "response-distr"){
                           plt <- ggplot2::ggplot(proj) +
                             ggplot2::geom_histogram(aes(x = !!sym(resp_var))) +
                             ggplot2::labs(x = resp_lab,
                                           y = "Count") +
                             ggplot2::theme_bw() +
                             basic_thm
                           
                           if ("Year" %in% colnames(proj)){
                             plt <- plt +
                               ggplot2::facet_wrap(~Year)
                           }
                         }
                         
                         # distribution of response x sd... - heatmap
                         if (type == "bivariate-distr"){
                           plt <- ggplot2::ggplot(proj) +
                             ggplot2::geom_hex(aes(x = !!sym(resp_var),
                                                   y = !!sym(sd_var))) +
                             ggplot2::labs(x = resp_lab,
                                           y = sd_lab,
                                           fill = "Count") +
                             ggplot2::scale_fill_gradient(low = "gainsboro",
                                                          high = "dodgerblue3") +
                             # ggplot2::scale_colour_gradient(low = "gainsboro",
                             #                                high = "dodgerblue3")
                             ggplot2::theme_bw() +
                             basic_thm
                           
                           if ("Year" %in% colnames(proj)){
                             plt <- plt +
                               ggplot2::facet_wrap(~Year)
                           }
                         }
                         
                         (plt)
                         
                       }
                       
                     ),
                     # Private entries to be created
                     private = list(
                       siteids = list(),
                       length = list()
                     )
)


