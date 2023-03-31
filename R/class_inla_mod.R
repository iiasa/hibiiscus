#' Template objects for inla model 
#' 
#' @description 
#' This is the template object for an inla model class.
#' It contains all the functions inherent to the model class. Any functions in 
#' the package that operate on or interact with the model(summarising, plottting) 
#' need to be coded in here or provided separately.
#' @keywords internal
#' @import R6
#' @noRd
hib_inla_model <- R6Class("hib_inla_model",
                          # Public data objects and functions
                          public = list(
                            version = NULL,
                            last_update = NULL,
                            # paths = NULL, # ??
                            data = NULL, # - restrict data to that used in model fitting... + version?
                            model = NULL,
                            # fixedef_coefficients = NULL,
                            # metrics = NULL,
                            processing_flags = NULL,
                            # TODO: what else to include
                            # 
                            # inla output is fairly basic, but extensive...
                            # need to have main effects
                            # important to store random effects params/hypers
                            user_formula = NULL,
                            inla_fitted_formula = NULL,
                            inla_incl_spatial_term = NULL,
                            inla_stack = NULL,
                            inla_mesh = NULL,
                            inla_mesh_settings = NULL,
                            
                            # Executed when created
                            # TODO
                            # Need to store both user provided (original), and  inla-input formula
                            # Need to store model stack (mesh etc...)
                            # need to store A matrix
                            # flag incl_spatial...
                            # store priors etc
                            initialize = function(db_obj, # ??
                                                  mod_df,
                                                  inla_obj,
                                                  # fixef_df,
                                                  # metr_vec,
                                                  # %%%
                                                  user_formula,
                                                  inla_fitted_formula,
                                                  inla_incl_spatial_term,
                                                  inla_stack,
                                                  inla_mesh,
                                                  inla_mesh_settings,
                                                  # %%%
                                                  n_spp
                                                  # set_up_specs
                                                  # spp df?
                            ){
                              # - some assertion checks
                              assertthat::assert_that(
                                # Check database obj
                                # db_obj$checkself(),
                                inherits(db_obj, "Database"),
                                # Check model data
                                is.data.frame(mod_df),
                                # Check model
                                class(inla_obj) == "inla",
                                # Check formulas
                                class(user_formula) == "formula",
                                class(inla_fitted_formula) == "formula",
                                # assertthat::assert_that(attr(user_formula, "class") == "formula"),
                                # assertthat::assert_that(attr(inla_fitted_formula, "class") == "formula"),
                                # Check inla specific data
                                inherits(inla_stack, "inla.data.stack"),
                                is.logical(inla_incl_spatial_term),
                                inherits(inla_mesh,'inla.mesh') || is.null(inla_mesh),
                                # Check coefficient df and metrics
                                # is.data.frame(fixef_df),
                                # is.numeric(metr_vec),
                                is.list(inla_mesh_settings),
                                is.numeric(n_spp)
                              )
                              
                              # ------ #
                              message(crayon::green("Initialising fitted hib model object ...") )
                              
                              # %%%%%%%%%%%%%%%
                              
                              # Save inputs in self
                              # Set mod_df as data
                              self$data <- mod_df
                              # Set model as inla_obj
                              self$model <- inla_obj
                              # Set fixed_coefficients as fixef_df
                              # self$fixedef_coefficients <- fixef_df
                              # Set metrics as metr_vec
                              # self$metrics <- metr_vec
                              
                              self$user_formula <- user_formula
                              self$inla_fitted_formula <- inla_fitted_formula
                              
                              self$inla_incl_spatial_term <- inla_incl_spatial_term
                              self$inla_stack <- inla_stack
                              self$inla_mesh <- inla_mesh
                              
                              # Incl processing flags from hib_db
                              self$processing_flags <- db_obj$processing_flags
                              
                              # Calc n sites/species used in  model - need to calc spp with div metrics!!
                              N_sites <- dplyr::n_distinct(mod_df$SSBS)
                              
                              self$processing_flags$N_spp <- n_spp
                              
                              self$version <- db_obj$version
                              self$last_update <- lubridate::today()
                              
                              # %%%%%%%%%%%%%%%
                              
                            },
                            # Functions
                            
                            # Base printout function
                            print = function(){
                              stats <- list(number_sites = dplyr::n_distinct(self$data$SSBS),
                                            # number_lu = length( unique(self$data$Predominant_land_use) ),
                                            # number_lui = length(unique(self$data$Use_intensity)),
                                            number_species = self$processing_flags$N_spp ) #,
                              # filters_applied = paste(self$processing_flags$filtering_steps, collapse = "; "),
                              # div_metrics = paste(sort(unique(self$diversity$Metric)), collapse = ", ") )  
                              
                              message(
                                paste0("Hibiiscus fitted inla model\n",
                                       "(Last updated: ", crayon::green( self$last_update ), ")",
                                       "\n",
                                       "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                                       # "\t Land use (class x intensity):  ", crayon::green(stats$number_lu), " x ",
                                       # crayon::green(stats$number_lui), "\n",
                                       "\t Number of species:             ", crayon::green(stats$number_species), "\n",
                                       "\n\n",
                                       "Summary of fitted model:\n"
                                )
                              )
                              print(summary(self$model))
                            },
                            # Save object
                            save_hib = function(output_path = NULL){
                              assertthat::assert_that(
                                is.null(output_path) || is.character(output_path)
                              )
                              
                              if(!assertthat::has_extension(output_path, "rds")){
                                output_path <- paste0(output_path, ".rds")
                              }
                              
                              if (file.exists(output_path)){
                                message(crayon::yellow(paste0("Warning: '", output_path, "' already exists.",
                                                              "\n",
                                                              "This file will be updated.\n")))
                              }
                              
                              message(crayon::green(paste0("Saving hib model to: '", output_path, "' ...\n")) )
                              
                              saveRDS(self, output_path)
                              assertthat::assert_that(file.exists(output_path),
                                                      msg = "Writing a copy of the model failed...")
                              invisible()
                            },
                            # TODO: plot function ...
                            plot = function(){
                              # sjPlot::plot_model(self$model, 
                              #                    show.intercept = T)
                              stop("Function in progress...")
                            }
                          ),
                          private = list()
)

