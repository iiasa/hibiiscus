#' Template objects for lme4 model 
#' 
#' @description 
#' This is the template object for an lme4 model class.
#' It contains all the functions inherent to the model class. Any functions in 
#' the package that operate on or interact with the model(summarising, plottting) 
#' need to be coded in here or provided separately.
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
                          inherits(db_obj, "Database"),
                          # Check model data
                          is.data.frame(mod_df),
                          # Check model
                          assertthat::assert_that(class(lme4_obj) == "lmerMod" || class(lme4_obj) == "glmerMod"), 
                          # Check coefficient df and metrics
                          is.data.frame(fixef_df),
                          is.numeric(metr_vec),
                          is.numeric(n_spp)
                          )
                        
                        # ------ #
                        message(crayon::green("Initialising fitted hib model object ...") )
                        
                        # %%%%%%%%%%%%%%%
                        
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
                        # N_sites <- dplyr::n_distinct(mod_df$SSBS)
                        
                        if ("SSBS" %in% colnames(mod_df)){
                          self$processing_flags$N_sites <- dplyr::n_distinct(mod_df$SSBS)
                        } else if ("SSBS1" %in% colnames(mod_df)){
                          # print("SSBS1 found")
                          self$processing_flags$N_sites <- length(unique(c(mod_df$SSBS1, mod_df$SSBS2)))
                        }
                        
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
                        stats <- list(number_sites = self$processing_flags$N_sites,
                                      # number_lu = length( unique(self$data$Predominant_land_use) ),
                                      # number_lui = length(unique(self$data$Use_intensity)),
                                      number_species = self$processing_flags$N_spp ) #,
                                      # filters_applied = paste(self$processing_flags$filtering_steps, collapse = "; "),
                                      # div_metrics = paste(sort(unique(self$diversity$Metric)), collapse = ", ") )  
                        # TODO:
                        # replicate INLA style printing!!
                        # summarise fixedeff coefs table
                        # display random effects types...
                        
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
                      plot = function(type = "fixed-effects"){
                        # sjPlot::plot_model(self$model, 
                        #                    show.intercept = T)
                        
                        assertthat::assert_that(is.character(type))
                        # stop("Function in progress...")
                        # Plot options
                        type = match.arg(type,
                                         c("fixed-effects",
                                           "random-effects"))
                        
                        basic_thm <- theme(axis.text = element_text(size = 16),
                                           axis.title = element_text(size = 20),
                                           plot.title = element_text(size = 20),
                                           plot.subtitle = element_text(size = 18),
                                           strip.text.x = element_text(size = 18),
                                           strip.text.y = element_text(size = 18),
                                           strip.background = element_blank(),
                                           legend.text = element_text(size = 16),
                                           legend.title = element_text(size = 18),
                                           legend.text.align = 0) 
                        
                        # 1. fixed effects - summary
                        if (type == "fixed-effects"){
                          fixef_df <- self$fixedef_coefficients
                          
                          # Set a baseline- Intercept
                          fixef_df$Coefficient_name <- factor(fixef_df$Coefficient_name)
                          if ("(Intercept)" %in% fixef_df$Coefficient_name){
                            fixef_df$Coefficient <- relevel(fixef_df$Coefficient_name,
                                                            ref = "(Intercept)")
                          }
                          
                          # %%%
                          plt <- ggplot(fixef_df) +
                            geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                            geom_errorbar(aes(x = Coefficient_name, 
                                              ymin = loCI,
                                              ymax = hiCI),
                                          width = 0) +
                            geom_point(aes(x = Coefficient_name,
                                           y = Coefficient_value)) +
                            
                            labs(y = "Estimate",
                                 x = element_blank(),
                                 title = "Fixed effect estimates") +
                            theme_bw() +
                            basic_thm +
                            coord_flip()
                          
                          plot(plt)
                          
                        }
                        #    fixed effects - posteriors
                        # 2. random iid
                        if (type == "random-effects"){
                          # if (is.null(self$user_formula$random)){
                          #   
                          # }
                          # Get random vars from random form
                          rand_vars <- all.vars(self$user_formula$random)
                          
                          rand_df <- as.data.frame(lme4::ranef(self$model)) %>%
                            dplyr::mutate(lwr = condval - 1.96*condsd, 
                                          upr = condval + 1.96*condsd) %>%
                            # order from small to large?
                            dplyr::group_by(grpvar) %>%
                            dplyr::arrange(condval, .by_group = T) %>%
                            dplyr::mutate(ID_ord = 1:n()) %>%
                            dplyr::ungroup() %>%
                            dplyr::mutate(sig_col = case_when(lwr > 0 ~ "+ve",
                                                              upr < 0 ~ "-ve",
                                                              T ~ "n.s."))
                          
                          # colour as +ve or -ve 
                          # plot - no names...
                          plt <- ggplot(rand_df) +
                            geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                            geom_errorbar(aes(x = ID_ord, 
                                              ymin = lwr,
                                              ymax = upr,
                                              colour = sig_col),
                                          show.legend = F,
                                          alpha = 0.5,
                                          width = 0) +
                            geom_point(aes(x = ID_ord,
                                           y = condval,
                                           colour = sig_col),
                                       show.legend = F,
                                       alpha = 0.7) +
                            labs(y = "Estimate",
                                 x = element_blank(),
                                 title = "Random effect estimates") +
                            scale_color_manual(values = c("+ve" = "dodgerblue3",
                                                          "-ve" = "firebrick3",
                                                          "n.s."= "grey30")) +
                            facet_wrap(~grpvar, scales = "free_x") +
                            theme_bw() +
                            basic_thm +
                            theme(axis.text.x = element_blank(),
                                  aspect.ratio = 1)
                          
                          plot(plt)
                        }
                      }
                    ),
                    private = list()
)