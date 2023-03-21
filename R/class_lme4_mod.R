#' Template objects for databases 
#' 
#' @description 
#' This is the template object for a database containing PREDICTS or PREDICTS-like data.
#' It contains all the functions inherent to the database. Any functions in the package that is 
#' operating or interacting with the data (preparation, linking, modelling) need to be coded in here 
#' or provided separately.
#' @keywords internal
#' @import R6
#' @noRd
hib_lme4_model <- R6Class("hib_lme4_model",
                    # Public data objects and functions
                    public = list(
                      version = NULL,
                      last_update = NULL,
                      # paths = NULL, # ??
                      data = NULL, # - restrict data to that used in model fitting... + version?
                      model = NULL,
                      fixedef_coefficients = NULL,
                      metrics = NULL,
                      processing_flags = NULL,
                      # TODO: what else to include
                      # 
                      # lme4 output is pretty comprehensive, this is just a wrapper...
                      # need optimal ranef struct, if found
                      # need model ranks, if dredged..
                      
                      # Executed when created
                      initialize = function(db_obj, # ??
                                            mod_df,
                                            lme4_obj,
                                            fixef_df,
                                            metr_vec,
                                            n_spp
                                            # set_up_specs
                                            # spp df?
                                            # ???? - should this do the model fitting
                                            # or should that be in hib_fit_lme4
                                            # and this just return the resultant obj
                                            # facilitating summary, printing...
                                            ){
                        # - some assertion checks
                        assertthat::assert_that(
                          # Check database obj
                          # db_obj$checkself(),
                          # Check model data
                          is.data.frame(mod_df),
                          # Check model
                          assertthat::assert_that(class(lme4_obj) == "lmerMod" || class(lme4_obj) == "glmerMod"), 
                          # Check coefficient df and metrics
                          is.data.frame(fixef_df),
                          is.numeric(metr_vec),
                          is.numeric(n_spp)
                          )
                        
                        # TODO
                        # ------ #
                        message(crayon::green("Initialising fitted hib model object ...") )
                        
                        # %%%%%%%%%%%%%%%
                        # self$checkself(sites, species)
                        
                        # Save inputs in self
                        # Set mod_df as data
                        self$data <- mod_df
                        # Set model as lme4_obj
                        self$model <- lme4_obj
                        # Set fixed_coefficients as fixef_df
                        self$fixedef_coefficients <- fixef_df
                        # Set metrics as metr_vec
                        self$metrics <- metr_vec
                        
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
                        # N spp/sites...
                        # print(summary(self$model))
                        stats <- list(number_sites = dplyr::n_distinct(self$data$SSBS),
                                      # number_lu = length( unique(self$data$Predominant_land_use) ),
                                      # number_lui = length(unique(self$data$Use_intensity)),
                                      number_species = self$processing_flags$N_spp ) #,
                                      # filters_applied = paste(self$processing_flags$filtering_steps, collapse = "; "),
                                      # div_metrics = paste(sort(unique(self$diversity$Metric)), collapse = ", ") )  
                        
                        
                        message(
                          paste0("Hibiiscus fitted lme4 model\n",
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
                      # TODO: plot function - just use sjplot?
                      plot = function(){

                        sjPlot::plot_model(self$model, 
                                           show.intercept = T)
                      }
                    ),
                    private = list()
)