#' Template object for cross-validation outputs
#' 
#' @description 
#' This is the template object for a cross-validation class.
#' It contains all the functions inherent to the cross-validation class. 
#' Any functions in the package that operate on or interact with the class (summarising, plottting) 
#' need to be coded in here or provided separately.
#' @keywords internal
#' @import R6
#' @noRd
hib_cv_res <- R6Class("hib_cv_res",
                          # Public data objects and functions
                          public = list(
                            # version = NULL,
                            last_update = NULL,
                            model_type = NULL,
                            original_model_obj = NULL, # All model info, except data, and functions
                            cv_results = NULL,
                            insample_results = NULL,
                            
                            # Executed when created
                            initialize = function(mod_obj,
                                                  cv_res,
                                                  insamp_res
                                                  
                            ){
                              # - some assertion checks
                              assertthat::assert_that(
                                # Check model
                                inherits(mod_obj, "hib_lme4_model"),
                                # Check cv outputs
                                is.list(cv_res),
                                is.list(insamp_res)
                              )
                              
                              # ------ #
                              message(crayon::green("Initialising cross-validation object ...") )
                              
                              # %%%%%%%%%%%%%%%
                              
                              # Save inputs in self
                             
                              core_items <- c("last_update", "version",
                                              "model", "processing_flags")
                              if (inherits(mod_obj, "hib_lme4_model")){ 
                                core_items <- c(core_items, c("metrics", 
                                                              "ci_type", 
                                                              "fixedef_coefficients"))
                                self$model_type <- "hib_lme4_model"
                                }
                              if (inherits(mod_obj, "hib_inla_model")){ 
                                core_items <- c(core_items, c("user_formula",
                                                              "inla_fitted_formula", 
                                                              "inla_incl_spatial_term",
                                                              "inla_mesh_settings"))
                                self$model_type <- "hib_inla_model"
                                }
                              if (inherits(mod_obj, "hib_brms_model")){ 
                                core_items <- c(core_items, c("metrics", 
                                                              "fixedef_coefficients"))
                                self$model_type <- "hib_inla_model"
                              }
                              
                              
                              self$original_model_obj <- purrr::map(.x = core_items, 
                                                                    .f = function(x, 
                                                                                  obj = mod_obj){
                                                                      obj[[x]]
                                                                    })
                              
                              names(self$original_model_obj) <- core_items
                              
                              self$cv_results <- cv_res
                              self$insample_results <- insamp_res
                              
                              self$last_update <- lubridate::today()
                              
                              # %%%%%%%%%%%%%%%
                              
                            },
                            # Functions
                            
                            # Base printout function
                            print = function(){
                              # Padding numbers function, from:
                              # https://github.com/jthomasmock/gtExtras/blob/HEAD/R/pad_fn.R
                              pad_fn <- function(x, nsmall = 2, pad0) {
                                # round and format values as text with specific number of decimals
                                round_x <- round(x, digits = nsmall)
                                fmt_x <- format(round_x, nsmall = nsmall)
                                
                                # calc number of trailing zeros
                                zero_len <- nchar(fmt_x) - nchar(sub("0*$", "", fmt_x))
                                
                                # create string of zero
                                if (pad0) {
                                  rep_zero <- strrep(0, zero_len)
                                } else {
                                  rep_zero <- ""
                                }
                                
                                # remove trailing zeros by position
                                fmt_out <- substr(fmt_x, 1, nchar(fmt_x) - zero_len)
                                
                                # add the non-breaking spaces to the formatted values
                                filled_out <- paste0(fmt_out, rep_zero)
                                
                                filled_out
                              }
                              
                              # stats <- list(number_sites = self$original_model_obj$processing_flags$N_sites,
                              #               number_species = paste(self$original_model_obj$processing_flags$N_spp, collapse = ", ") ) #,
                              
                              
                              
                              # Model form
                              m_form <- format(stats::formula(self$original_model_obj$model))
                              
                              
                              # Cv types and n folds per type
                              # make into a df
                              # cv_type_info <- purrr::map_dfr(.x = self$cv_results,
                              #                 .f = function(x){
                              #                   fld_typ <- x$fold_type
                              #                   if (fld_typ == "Custom"){
                              #                     fld_typ <- paste0(fld_typ, " (", x$fold_col, ")")
                              #                   }
                              #                   n_fld <- x$n_folds
                              #                   
                              #                   return(data.frame(Fold_type = fld_typ,
                              #                                     N_folds   = n_fld))
                              #                 })
                              
                            
                              
                              message(
                                paste0("Hibiiscus cross-validated lme4 model\n",
                                       "(Last updated: ", crayon::green( self$last_update ), ")",
                                       # "\n",
                                       # "\t Cross-validation types:        ", crayon::green(stats$number_sites), "\n"
                                       # "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                                       # # "\t Land use (class x intensity):  ", crayon::green(stats$number_lu), " x ",
                                       # # crayon::green(stats$number_lui), "\n",
                                       # "\t Number of species:             ", crayon::green(stats$number_species), "\n",
                                       "\n\n",
                                       "--------------------\n\n",
                                       "Summary of cross-validated model:\n\n",
                                       "Formula:                  ", crayon::green( m_form ) , "\n"
                                )
                              )
                              
                              # In-sample...
                              message(
                                paste0("--------------------\n",
                                       "Cross-validation type:    ", crayon::green("In-sample"), "\n", 
                                       "No. of folds:             ", crayon::green("1"), "\n", 
                                       "Predictive R^2:           ", crayon::green(self$insample_results$insample_metrics["Predictive_R_squared"]) , "\n")
                              )
                              
                              # per cv type
                              for (x in self$cv_results){
                                fld_typ <- x$fold_type
                                if (fld_typ == "Custom"){
                                  fld_typ <- paste0(fld_typ, " (", x$fold_col, ")")
                                }
                                n_fld <- x$n_folds
                                                  
                                # corsq <- x$cv_metrics["Correlation_squared"]
                                # rmse <- x$cv_metrics["Root_mean_squared_error"]
                                r2 <- x$cv_metrics["Predictive_R_squared"]
                                
                                
                                message(
                                  paste0("--------------------\n",
                                         "Cross-validation type:    ", crayon::green(fld_typ), "\n", 
                                         "No. of folds:             ", crayon::green(n_fld), "\n", 
                                         "Predictive R^2:           ", crayon::green(r2) , "\n")
                                )
                              }
                              # print(cv_type_info)
                              
                              # maybe some original model metrics..
                              
                              # metr <- pad_fn(self$metrics, 3, F)
                              # names(metr) <- names(self$metrics)
                              
                              # message(paste0("\n\n",
                              #                # "--------------------\n",
                              #                # "Random effects:\n",
                              #                # ranef_chr, "\n\n",
                              #                "Performance metrics:\n",
                              #                "\tR-squared (marg.) \t ", metr["R2m"], "\n",
                              #                "\tR-squared (cond.) \t ", metr["R2c"], "\n",
                              #                "\tAIC               \t ", metr["AIC"], "\n",
                              #                "\tAICc              \t ", metr["AICc"], "\n",
                              #                "\tBIC               \t ", metr["BIC"], "\n"))
                              
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
                              
                              message(crayon::green(paste0("Saving cross-validation results to: '", output_path, "' ...\n")) )
                              
                              saveRDS(self, output_path)
                              assertthat::assert_that(file.exists(output_path),
                                                      msg = "Writing a copy of the results failed...")
                              invisible()
                            },
                            
                            # TODO: plot function - 
                            plot = function(){
                              # type = "coefficients"
                              stop("Plotting functionality still in development...")
                              # Types:
                              #   -coefs - insample v cv (all folds v summary)
                              #   -metrics - over all fold types ; 
                              #              if multiple metrics, facter_wrap
                              #   -obs v pred, facet by fold_type, colour by fold_id
                              
                              assertthat::assert_that(is.character(type))
                              # Plot options
                              # type = match.arg(type,
                              #                  c("",
                              #                    ""))
                              
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
                              
                              
                              # (plt)
                            }
                          ),
                          private = list()
)
