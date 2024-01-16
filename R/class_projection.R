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
HibProjection <- R6Class("HibProjection",
  # TODO: @export
                              # hib_mod_projection "HibProjection"
                              # Public data objects and functions
                              public = list(

                                #' @field last_update (`character`)\cr
                                #' Identifies the most recent update to the object (YYYY-MM-DD).
                                last_update = NULL,
                                
                                #' @field model_type (`character`)\cr
                                #' Identifies the type of model; one of "HibModelLME4", 
                                #' "HibModelINLA", "HibModelBRMS".
                                model = NULL, # should this be a whole model class??
                                
                                #' @param projection (`data.frame`, `SpatRaster`, or `list`)\cr
                                #' The results of model projection.
                                projection = NULL,
                                
                                #' @param projection_specs (`list`)\cr
                                #' List of specifications used for modelling.
                                projection_specs = NULL,
                                
                                #' @param aoa (`list`)\cr
                                #' Results of Area-of-Applicability calculations.
                                aoa = NULL,
                                
                                # processing_flags = NULL,
                                
                                #' @description
                                #' Creates a new instance of this [R6][R6::R6Class] class.
                                #' Note, this object is constructed via the [hib_project_model]
                                #' function found within this package.
                                #' @param obj (`HibModelXXXX`)\cr
                                #'   The model object used in cross-validation.
                                #' @param proj (`data.frame`, `SpatRaster`, or `list`)\cr
                                #'   The results of model projection.
                                #' @param rescale_pred (`logical`)\cr
                                #'   Describes whether projected values are rescaled
                                #'   relative to the baseline land-use class.
                                #' @param lme4_specs (`list`)\cr
                                #'   List of specifications used for lme4 models.
                                #' @param inla_specs (`list`)\cr
                                #'   List of specifications used for inla models.
                                #' @param brms_specs (`list`)\cr
                                #'   List of specifications used for brms models.
                                #' @param aoa_res (`list`, or `NULL`)\cr
                                #'   List of Area-of-Applicability results.
                                # Executed when created
                                initialize = function(obj,
                                                      proj,
                                                      backtransform,
                                                      rescale_pred,
                                                      lme4_specs,
                                                      inla_specs,
                                                      brms_specs,
                                                      aoa_res){
                                  # Check that input is either a character filepath
                                  # Or dataframe/tibble
                                  assertthat::assert_that(
                                    inherits(obj, "HibModelLME4") || inherits(obj, "HibModelINLA") || inherits(obj, "HibModelBRMS"), # obj is a hib model
                                    is.data.frame(proj) || inherits(proj, "SpatRaster") || (is.list(proj) & !is.data.frame(proj)), # or a list of SpatRasters
                                    is.function(backtransform) || is.null(backtransform),
                                    is.logical(rescale_pred),
                                    is.list(lme4_specs) || is.null(lme4_specs),
                                    is.list(inla_specs) || is.null(inla_specs),
                                    is.list(brms_specs) || is.null(brms_specs),
                                    is.list(aoa_res) || is.null(aoa_res)
                                  )
                                  
                                  # ------ #
                                  message(crayon::green("\nInitialising HibProjection object ...") )
                                  
                                  # Store projection object (list of data)
                                  self$last_update <- lubridate::now()
                                  
                                  self$model <- obj
                                  self$projection <- proj
                                  proj_specs <- list("Backransformation_function" = backtransform,
                                                     "Rescaled_prediction" = rescale_pred)
                                  
                                  if (inherits(obj, "HibModelLME4")){
                                    self$projection_specs <- c(proj_specs, lme4_specs)
                                  }
                                  if (inherits(obj, "HibModelINLA")){
                                    # Package check
                                    if (!requireNamespace("INLA", quietly = TRUE)) {
                                      stop(
                                        "Package \"INLA\" must be installed to use this function.",
                                        call. = FALSE
                                      )
                                    }
                                    # OR
                                    # rlang::check_installed("INLA", reason = "to use INLA-based model fitting.")
    
                                    self$projection_specs <- c(proj_specs, inla_specs)
                                  }
                                  if (inherits(obj, "HibModelBRMS")){
                                    self$projection_specs <- c(proj_specs, brms_specs)
                                  }
                                  
                                  self$aoa <- aoa_res
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
                                
                                #' @description
                                #' Base printout function presents core information about the 
                                #' object.
                                #'
                                #' @return (`NULL`).
                                print = function() {
                                  # stop("Print functionality in progress...")
                                  # Compute some stats/get some info
                                  
                                  mod_type <- ""
                                  pred_vars <- ""
                                  n_pred_vars <- NA
                                  if (inherits(self$model, "HibModelLME4")){
                                    mod_type <- "lme4"
                                    
                                    all_covars <- all.vars(formula(delete.response(terms(formula(self$model$model)))))
                                    rand_vars <- names(lme4::ranef(self$model$model))
                                    fix_vars <- all_covars[!all_covars %in% rand_vars]
                                  }
                                  
                                  if (inherits(self$model, "HibModelINLA")){
                                    mod_type <- "INLA"
                                    
                                    fix_vars <- all.vars(formula(delete.response(terms(self$model$user_formula$response_fixed))))
                                    fix_vars <- fix_vars[!fix_vars == "Intercept"]
                                    fix_vars <- fix_vars[!fix_vars == "(Intercept)"]
                                    
                                    if (self$model$inla_incl_spatial_term){
                                      fix_vars <- c(fix_vars, "spatial-random-effect")
                                    }
                                  }
                                  
                                  if (inherits(self$model, "HibModelBRMS")){
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
                                    paste0("Hibiiscus projection\n",
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
                                
                                # #' @description
                                # #' Function to save the object to the specified file-path 
                                # #' ('output_path', must be '.rds').
                                # #'
                                # #' @param output_path (`character()`)
                                # #' @return (`NULL`).
                                # # TODO - have a generic save_hib function for all hib datatypes
                                # # Everything is just an RDS anyway
                                # save_projection = function(output_path = NULL){
                                #   # TODO: check this works with a terra-based projection
                                #   # terra objects are pointer-associated...
                                #   assertthat::assert_that(
                                #     is.character(output_path)
                                #   )
                                  
                                #   if(!assertthat::has_extension(output_path, "rds")){
                                #     output_path <- paste0(output_path, ".rds")
                                #   }
                                  
                                #   if (!dir.exists(dirname(output_path))){
                                #     message(crayon::yellow(paste0("\nWarning: '", dirname(output_path), "' doesn't exist.",
                                #                                   "\n",
                                #                                   "This folder will be created.\n")))
                                #     dir.create(dirname(output_path), recursive = T)
                                #   }
                                  
                                #   if (file.exists(output_path)){
                                #     message(crayon::yellow(paste0("\nWarning: '", output_path, "' already exists.",
                                #                                   "\n",
                                #                                   "This file will be updated.\n")))
                                #   }
                                  
                                #   message(crayon::green(paste0("Saving projection to: '", output_path, "'\n")) )
                                  
                                #   saveRDS(self, output_path)
                                #   assertthat::assert_that(file.exists(output_path),
                                #                           msg = "Writing a copy of the projection failed...")
                                #   invisible()
                                # },
                                
                                #' @description
                                #' Function to generate figures summarising the
                                #' projection results.
                                #' 
                                #' @param type (`character()`) Options are: 
                                #' "response-map", "uncertainty-map", "bivariate-map",
                                #' "response-distr", "uncertainty-distr", "bivariate-distr".
                                #' @param response_var (`character`, or `NULL`) 
                                #' @param uncertainty_var (`character`, or `NULL`) 
                                #' @param aoa_mask (`logical`) 
                                #' @param aoa_var (`character`, or `NULL`) 
                                #' @return (`ggplot`).
                                plot = function(type = "response-map",
                                                response_var = NULL,
                                                uncertainty_var = NULL,
                                                aoa_mask = FALSE,
                                                aoa_var = NULL
                                                # aoa_clr ??
                                                ){
                                  
                                  assertthat::assert_that(
                                    is.character(type),
                                    is.character(response_var) || is.null(response_var),
                                    is.character(uncertainty_var) || is.null(uncertainty_var),
                                    is.logical(aoa_mask),
                                    is.character(aoa_var) || is.null(aoa_var)
                                  )
                                  
                                  # add response_var, and error_var (or uncert_var) - for specifying plot details
                                  type <- match.arg(type,
                                                    c("response-map", "uncertainty-map", "bivariate-map",
                                                      "response-distr", "uncertainty-distr", "bivariate-distr"))
                                  
                                  
                                  # TODO: incorporate aoa info
                                  # Colour "-distr" plots - or rm from plot??
                                  
                                  
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
                                  
                                  # Get projection info - self$projection
                                  proj <- self$projection
                                  if (aoa_mask){
                                    aoa_exists <- !is.null(self$aoa)
                                    if (!aoa_exists){
                                      message(crayon::yellow(paste0("\nWarning: AOA result does not exist, ",
                                                                    "AOA cannot be plotted.")))
                                      aoa_mask <- F
                                      aoa_var <- NULL
                                    }

                                    if (aoa_exists){
                                      aoa <- self$aoa$proj_aoa  
                                    }
                                  }
                                  

                                  # TODO - incl respnse var/error var in args
                                  # get these from data accordingly...
                                  
                                  # Possible response vars
                                  # (rescaled_)prediction_mean, (rescaled_)prediction_median,
                                  # (rescaled_)prediction_q2.5, (rescaled_)prediction_q97.5
                                  # INLA only: (rescaled_)prediction_mode
                                  
                                  # Possible error/uncert vars
                                  # (rescaled_)prediction_sd, 
                                  # BRMS only: (rescaled_)prediction_mad
                                  
                                  # neaten these for axis/plot labels - 
                                  
                                  
                                  # INLA - (rescaled_)prediction_sd
                                  # BRMS - (rescaled_)prediction_error
                                  
                                  # Note: if proj is a list of rasters, assertthat::hasname doesn't work
                            
                                  
                                  # Check specified response and uncert vars are in projection data
                                  # Do this after conversion to df...
                                  # if ((!is.null(response_var) & type %in% c("response-map", 
                                  #                                           "bivariate-map",
                                  #                                           "response-distr", 
                                  #                                           "bivariate-distr")) ){
                                  #   assertthat::assert_that(
                                  #     assertthat::has_name(proj, response_var),
                                  #     msg = "Error: specified response variable is not available in projection output."
                                  #   )
                                  # }
                                  # if ((!is.null(uncertainty_var) & type %in% c("uncertainty-map", 
                                  #                                              "bivariate-map",
                                  #                                              "uncertainty-distr", 
                                  #                                              "bivariate-distr")) ){
                                  #   assertthat::assert_that(
                                  #     assertthat::has_name(proj, uncertainty_var),
                                  #     msg = "Error: specified uncertainty variable is not available in projection output."
                                  #   )
                                  # }
                                  
                                  # If response/uncert are not provided, default to prediction_mean and prediction_sd
                                  if ((is.null(response_var) & type %in% c("response-map", 
                                                                           "bivariate-map",
                                                                           "response-distr", 
                                                                           "bivariate-distr")) ){
                                    message(crayon::yellow(paste0("\nWarning: 'response_var' not specified, ",
                                                                  "using 'prediction_mean' as a default.")))
                                    response_var <- "prediction_mean"
                                    # assertthat::assert_that(
                                    #   assertthat::has_name(proj, response_var),
                                    #   msg = "Error: specified response variable is not available in projection output."
                                    # )
                                  }
                                  if ((is.null(uncertainty_var) & type %in% c("uncertainty-map", 
                                                                              "bivariate-map",
                                                                              "uncertainty-distr", 
                                                                              "bivariate-distr")) ){
                                    message(crayon::yellow(paste0("\nWarning: 'uncertainty_var' not specified, ",
                                                                  "using 'prediction_sd' as a default.")))
                                    uncertainty_var <- "prediction_sd"
                                    # assertthat::assert_that(
                                    #   assertthat::has_name(proj, uncertainty_var),
                                    #   msg = "Error: specified uncertainty variable is not available in projection output."
                                    # )
                                  }
                                  
                                  # Set labels - tidy versions of response/uncertainty_var
                                  
                                  if (!is.null(response_var)){
                                    resp_lab <- gsub("_", " ", response_var) |>
                                      tools::toTitleCase()
                                    
                                    # replace "Mean" with "(mean)"
                                    resp_lab <- gsub("Mean", "(Mean)", resp_lab)
                                    # replace "Median" with "(median)"
                                    resp_lab <- gsub("Median", "(Median)", resp_lab)
                                    # replace "Q2.5" with "(2.5% quantile)"
                                    resp_lab <- gsub("Q2.5", "(2.5% quantile)", resp_lab)
                                    # replace "Q97.5" with "(97.5% quantile)"
                                    resp_lab <- gsub("Q97.5", "(97.5% quantile)", resp_lab)
                                    
                                    # TODO:
                                  }
                                  
                                  
                                  if (!is.null(uncertainty_var)){
                                    uncert_lab <- gsub("_", " ", uncertainty_var) |>
                                      tools::toTitleCase()
                                    
                                    # replace "Sd", with "(SD)"
                                    uncert_lab <- gsub("Sd", "SD", uncert_lab)
                                    # replace "Mad", with "(MAD)"
                                    uncert_lab <- gsub("Mad", "MAD", uncert_lab)
                                    
                                    # TODO:
                                    if (uncertainty_var == "prediction_sd"){
                                      uncert_lab <- "Uncertainty (SD)"
                                    }
                                    
                                    if (uncertainty_var == "prediction_mad"){
                                      uncert_lab <- "Uncertainty (MAD)"
                                    }
                                  }
                                  
                                  
                                  # TODO: check that resp and uncert vars are in proj
                                  # If lme4 without boot, no uncert...
                                  
                                  vars_ <- c(response_var, uncertainty_var)
                                  
                                  # If not rescaled..., just use a basic colour gradient
                                  # Initially, basic ggplot
                                  
                                  # TODO: add other plotting options...
                                  # such as...
                                  
                                  if (is.data.frame(proj) ){
                                    # terra::plot(proj$`2010`[[resp_var]])
                                    # TODO
                                    if (!all(vars_ %in% colnames(proj))){
                                      missing_vars <- vars_[!vars_ %in% colnames(proj)]
                                      message(crayon::red(paste("\nVariable/s missing:",
                                                                paste(missing_vars, collapse = " ; "),
                                                                "\n")))
                                      stop("Variable/s required for plotting are missing, see above for details.",
                                           call. = F)
                                    }
                                    
                                    vars__ <- c("x", "y", "Year", vars_)
                                    vars__ <- vars__[vars__ %in% colnames(proj)]
                                    proj <- proj[vars__]
                                  }
                                  if (inherits(proj, "SpatRaster") ){
                                    # terra::plot(proj$`2010`[[resp_var]])
                                    if (!all(vars_ %in% names(proj))){
                                      missing_vars <- vars_[!vars_ %in% names(proj)]
                                      message(crayon::red(paste("\nVariable/s missing:",
                                                                paste(missing_vars, collapse = " ; "),
                                                                "\n")))
                                      stop("Variable/s required for plotting are missing, see above for details.",
                                           call. = F)
                                    }
                                    
                                    proj <- terra::as.data.frame(proj[[vars_]], xy = T)

                                    if (aoa_mask){
                                      aoa <- terra::as.data.frame(aoa, xy = T)  
                                    }
                                    
                                  }
                                  if ((is.list(proj) & !is.data.frame(proj))){
                                    proj <- purrr::map(proj,
                                                       function(z){
                                                         if (!all(vars_ %in% names(z))){
                                                           missing_vars <- vars_[!vars_ %in% names(z)]
                                                           message(crayon::red(paste("\nVariable/s missing:",
                                                                                     paste(missing_vars, collapse = " ; "),
                                                                                     "\n")))
                                                           stop("Variable/s required for plotting are missing, see above for details.",
                                                                call. = F)
                                                         }
                                                         
                                                         terra::as.data.frame(z[[vars_]], xy = T)}) |>
                                      dplyr::bind_rows(.id = "Year")

                                    if (aoa_mask){
                                      aoa <- purrr::map(aoa,
                                                         function(z){
                                                          terra::as.data.frame(z, xy = T)
                                                          }) |>
                                        dplyr::bind_rows(.id = "Year")
                                    }
                                  }
                                  
                                  
                                  # TODO:
                                  # check that aoa_var is in colnames(proj_df)
                                  if (aoa_mask & !is.null(aoa_var)){
                                    assertthat::assert_that(
                                      assertthat::has_name(aoa, aoa_var)
                                    )
                                  }
                                  if (aoa_mask & is.null(aoa_var)){
                                    # get a vector of possible aoa cols
                                    aoa_cols <- colnames(aoa)
                                    aoa_cols <- aoa_cols[grepl("AOA_", aoa_cols)]
                                    assertthat::assert_that(
                                      length(aoa_cols) > 0
                                    )
                                    # randomly use 1...
                                    aoa_var <- aoa_cols[1]
                                    message(crayon::yellow(paste0("\nWarning: 'aoa_var' not specified, ",
                                                                  "using '", 
                                                                  aoa_var,
                                                                  "'.")))
                                  }


                                  # TODO: make use of tidyterra??
                                  # print(head(proj))
                                  
                                  
                                  if (type == "response-map"){
                                    # TODO: check x/y exist

                                    plt <- ggplot2::ggplot(proj) +
                                      ggplot2::geom_tile(ggplot2::aes(x = x, y = y,
                                                                      fill = !!dplyr::sym(response_var),
                                                                      colour = !!dplyr::sym(response_var))) +
                                      ggplot2::labs(fill = resp_lab,
                                                    colour = resp_lab) +
                                      ggplot2::coord_fixed() +
                                      ggplot2::theme_bw() +
                                      basic_thm
                                    
                                    # If rescaled, can also improve colour scheme 
                                    # 1 as midpoint, lower "bad" (0), higher "better" (max)
                                    if (grepl("rescaled_", response_var, ignore.case = T)){
                                      # if (self$projection_specs$Rescaled_prediction){
                                      max_resp <- max(proj[,response_var])
                                      
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
                                      
                                    } else {
                                      # if (self$projection_specs$Rescaled_prediction == F)
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
                                    # TODO: check x/y exist
                                    # TODO: add colour scale options...
                                    # c("Bluegill", "BlueGold", "BlueOr", "BlueYl", "Brown",
                                    #   "Brown2", "DkBlue","DkBlue2", "DkCyan","DkCyan2", "DkViolet",
                                    #   "DkViolet2", "GrPink","GrPink2", "PinkGrn", "PurpleGrn", "PurpleOr")
                                    legend <- biscale::bi_legend(pal = "GrPink2",
                                                                 dim = 3,
                                                                 xlab = paste("Larger", resp_lab),
                                                                 ylab = paste("Larger", uncert_lab),
                                                                 size = 12)
                                    
                                    # Normalise response and sd...??
                                    
                                    proj$var1 <- proj[[response_var]]
                                    proj$var2 <- proj[[uncertainty_var]]
                                    
                                    suppressWarnings(
                                      proj <- biscale::bi_class(proj, x = var1, y = var2, dim = 3, style = "quantile")
                                    )
                                    
                                    plt <- ggplot2::ggplot() +
                                      ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = bi_class),
                                                         data = proj,
                                                         show.legend = F) +
                                      biscale::bi_theme(base_size = 16) +
                                      biscale::bi_scale_fill(pal = "GrPink2", na.value = "transparent") +
                                      # ggplot2::labs(
                                        # title = ggplot2::element_blank() resp_lab
                                        # x = ggplot2::element_blank(),
                                        # y = ggplot2::element_blank()
                                      # ) +
                                      ggplot2::theme_bw() +
                                      basic_thm +
                                      ggplot2::theme(#legend.position = "none",
                                                     aspect.ratio = 1)
                                    
                                    if ("Year" %in% colnames(proj)){
                                      plt <- plt +
                                        ggplot2::facet_wrap(~Year)
                                    }
                                    
                                    # print(grepl("-map", type) && aoa_mask)
                                    if (aoa_mask){
                                      aoa <- aoa |>
                                        dplyr::filter(!!dplyr::sym(aoa_var) == 0)
                                      plt <- plt + 
                                        ggplot2::geom_tile(ggplot2::aes(x = x, y = y,
                                                                        alpha = "1"),
                                                          fill = "#FF00FF",# "purple",#"grey10",
                                                          colour = "#FF00FF",
                                                          data = aoa) +
                                        ggplot2::scale_alpha_manual(values = c("1" = 1),
                                                          labels = c("1" = gsub("AOA_", "CV type: ", aoa_var))) +
                                        ggplot2::labs(alpha = "Outside AOA") +
                                        ggplot2::theme(legend.position = "bottom")
                                    }

                                    plt <- plt +
                                      ggplot2::labs(title = "Model projection: Bivariate map")
                                    

                                    plt <- cowplot::ggdraw() +
                                      cowplot::draw_plot(plt, 0, 0, 0.75, 1) +
                                      cowplot::draw_plot(legend, 0.75, .4, 0.25, 0.25)
                                  }
                                  
                                  # "uncertainty-map"
                                  if (type == "uncertainty-map"){
                                    # TODO: check x/y exist
                                    # %%%
                                    plt <- ggplot2::ggplot(proj) +
                                      ggplot2::geom_tile(ggplot2::aes(x = x, y = y,
                                                                      fill = !!dplyr::sym(uncertainty_var),
                                                                      colour = !!dplyr::sym(uncertainty_var))) +
                                      ggplot2::labs(fill = uncert_lab,
                                                    colour = uncert_lab) +
                                      ggplot2::coord_fixed() +
                                      ggplot2::theme_bw() +
                                      basic_thm
                                    
                                    plt <- plt +
                                      ggplot2::scale_fill_gradient(low = "gainsboro",
                                                                   high = "dodgerblue3") +
                                      ggplot2::scale_colour_gradient(low = "gainsboro",
                                                                     high = "dodgerblue3")
                                    
                                    
                                    # If rescaled, can also improve colour scheme 
                                    # 1 as midpoint, lower "bad" (0), higher "better" (max)
                                    # if (grepl("Rescaled_", response_var)){
                                    #   # if (self$projection_specs$Rescaled_prediction){
                                    #   max_resp <- max(proj[,response_var])
                                    #   
                                    #   plt <- plt +
                                    #     ggplot2::scale_fill_gradient2(low = "firebrick3",
                                    #                                   mid = "gainsboro",
                                    #                                   high = "dodgerblue3",
                                    #                                   midpoint = 1,
                                    #                                   limits = c(0,max_resp)) +
                                    #     ggplot2::scale_colour_gradient2(low = "firebrick3",
                                    #                                     mid = "gainsboro",
                                    #                                     high = "dodgerblue3",
                                    #                                     midpoint = 1,
                                    #                                     limits = c(0,max_resp))
                                    #   
                                    # } else {
                                    #   # if (self$projection_specs$Rescaled_prediction == F)
                                    #   plt <- plt +
                                    #     ggplot2::scale_fill_gradient(low = "gainsboro",
                                    #                                  high = "dodgerblue3") +
                                    #     ggplot2::scale_colour_gradient(low = "gainsboro",
                                    #                                    high = "dodgerblue3")
                                    # }
                                    
                                    if ("Year" %in% colnames(proj)){
                                      plt <- plt +
                                        ggplot2::facet_wrap(~Year)
                                    }
                                    # %%%
                                    
                                  }

                                  if (type %in% c("response-map", "uncertainty-map") & 
                                      aoa_mask){
                                    aoa <- aoa |>
                                      dplyr::filter(!!dplyr::sym(aoa_var) == 0)
                                    plt <- plt + 
                                      ggplot2::geom_tile(ggplot2::aes(x = x, y = y,
                                                                      alpha = "1"),
                                                        fill = "#FF00FF",
                                                        colour = "#FF00FF",
                                                        data = aoa) +
                                      ggplot2::scale_alpha_manual(values = c("1" = 1),
                                                        labels = c("1" = gsub("AOA_", "CV type: ", aoa_var))) +
                                      ggplot2::labs(alpha = "Outside AOA") +
                                      ggplot2::guides(alpha = ggplot2::guide_legend(order = 2),
                                                      fill = ggplot2::guide_colourbar(order = 1),
                                                      colour = ggplot2::guide_colourbar(order = 1))
                                  }
                                  
                                  # distribution of response... - histogram
                                  if (type == "response-distr"){
                                    # TODO: replicate this in each -distr plt
                                    if (aoa_mask){
                                      # TODO: check x/y exist
                                      df <- proj |>
                                        dplyr::left_join(aoa,
                                                          by = c("x", "y"))
                                      plt <- ggplot2::ggplot(df) +
                                        ggplot2::geom_histogram(ggplot2::aes(x = !!dplyr::sym(response_var), 
                                                                             fill = as.character(!!dplyr::sym(aoa_var)) ),
                                                                alpha = 0.7) +
                                        ggplot2::scale_fill_manual(values = c("0" = "#FF00FF",
                                                                     "1" = "grey35"),
                                                          labels = c("0" = "Outside",
                                                                     "1" = "Inside")) +
                                        ggplot2::labs(x = resp_lab,
                                                    y = "Count",
                                                    fill = paste0("AOA ",
                                                                  gsub("AOA_", "\nCV type: ", aoa_var)) ) +
                                      ggplot2::theme_bw() +
                                      basic_thm

                                    } else { 
                                      plt <- ggplot2::ggplot(proj) +
                                        ggplot2::geom_histogram(ggplot2::aes(x = !!dplyr::sym(response_var)),
                                                                alpha = 0.7) +
                                        ggplot2::labs(x = resp_lab,
                                                      y = "Count") +
                                        ggplot2::theme_bw() +
                                        basic_thm
                                    }
                                    
                                    if ("Year" %in% colnames(proj)){
                                      plt <- plt +
                                        ggplot2::facet_wrap(~Year)
                                    }
                                  }
                                  
                                  # distribution of response x sd... - heatmap
                                  if (type == "bivariate-distr"){
                                    # if (aoa_mask){...}
                                    # AOA - maybe have side-by side plots? - facet by inside/outside AOA
                                    if (aoa_mask){
                                      # TODO: check x/y exist
                                      df <- proj |>
                                        dplyr::left_join(aoa,
                                                          by = c("x", "y")) |>
                                        dplyr::mutate(AOA = !!dplyr::sym(aoa_var))

                                      plt <- ggplot2::ggplot(df) +
                                        ggplot2::geom_hex(ggplot2::aes(x = !!dplyr::sym(response_var),
                                                                       y = !!dplyr::sym(uncertainty_var))) +
                                                                       # colour = as.character(AOA))) +
                                        ggplot2::labs(x = resp_lab,
                                                      y = uncert_lab,
                                                      fill = "Count") +
                                        ggplot2::scale_fill_gradient(low = "gainsboro",
                                                                     high = "dodgerblue3") +
                                        # ggplot2::scale_colour_manual(values = c("0" = "#FF00FF",
                                        #                                         "1" = NA)) +
                                        # ggplot2::scale_colour_gradient(low = "gainsboro",
                                        #                                high = "dodgerblue3")
                                        ggplot2::facet_wrap(~AOA,
                                                            labeller = ggplot2::labeller("AOA" = c("0" = "Outside AOA",
                                                                                                    "1" = "Inside AOA"))) +
                                        ggplot2::theme_bw() +
                                        basic_thm

                                      if ("Year" %in% colnames(df)){
                                        plt <- plt +
                                          ggplot2::facet_wrap(AOA~Year,
                                                              labeller = ggplot2::labeller("AOA" = c("0" = "Outside AOA",
                                                                                                    "1" = "Inside AOA")))
                                      }

                                    } else {
                                      plt <- ggplot2::ggplot(proj) +
                                        ggplot2::geom_hex(ggplot2::aes(x = !!dplyr::sym(response_var),
                                                                       y = !!dplyr::sym(uncertainty_var))) +
                                        ggplot2::labs(x = resp_lab,
                                                      y = uncert_lab,
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
                                  }

                                  # "uncertainty-distr",
                                  if (type == "uncertainty-distr"){
                                    # TODO: replicate this in each -distr plt
                                    if (aoa_mask){
                                      # TODO: check x/y exist
                                      df <- proj |>
                                        dplyr::left_join(aoa,
                                                          by = c("x", "y"))
                                      plt <- ggplot2::ggplot(df) +
                                        ggplot2::geom_histogram(ggplot2::aes(x = !!dplyr::sym(uncertainty_var), 
                                                                             fill = as.character(!!dplyr::sym(aoa_var)) ),
                                                                alpha = 0.7) +
                                        ggplot2::scale_fill_manual(values = c("0" = "#FF00FF",
                                                                     "1" = "grey35"),
                                                          labels = c("0" = "Outside",
                                                                     "1" = "Inside")) +
                                        ggplot2::labs(x = uncert_lab,
                                                    y = "Count",
                                                    fill = paste0("AOA ",
                                                                  gsub("AOA_", "\nCV type: ", aoa_var)) ) +
                                      ggplot2::theme_bw() +
                                      basic_thm

                                    } else { 
                                      plt <- ggplot2::ggplot(proj) +
                                        ggplot2::geom_histogram(ggplot2::aes(x = !!dplyr::sym(uncertainty_var)),
                                                                alpha = 0.7) +
                                        ggplot2::labs(x = uncert_lab,
                                                      y = "Count") +
                                        ggplot2::theme_bw() +
                                        basic_thm
                                    }
                                    # plt <- ggplot2::ggplot(proj) +
                                    #   ggplot2::geom_histogram(ggplot2::aes(x = !!dplyr::sym(uncertainty_var))) +
                                    #   ggplot2::labs(x = uncert_lab,
                                    #                 y = "Count") +
                                    #   ggplot2::theme_bw() +
                                    #   basic_thm
                                    
                                    if ("Year" %in% colnames(proj)){
                                      plt <- plt +
                                        ggplot2::facet_wrap(~Year)
                                    }
                                  }

                                  
                                  # Convert type to a title
                                  type_title <- gsub("-", "_", type)
                                  type_title <- tools::toTitleCase(type_title)
                                  type_title <- gsub("_", " ", type_title)
                                  type_title <- gsub("distr", "distribution", type_title)
                                  type_title <- paste0("Model projection: ", type_title)

                                  plt <- plt +
                                    ggplot2::labs(title = type_title)
                                  
                                  (plt)
                                  
                                }
                                
                              ),
                              # Private entries to be created
                              private = list(
                                siteids = list(),
                                length = list()
                              )
)

