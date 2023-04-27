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
                                class(user_formula) == "formula" || is.list(user_formula),
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
                              # inla_obj$summary.fixed
                              # Set metrics as metr_vec
                              # self$metrics <- metr_vec
                              
                              self$user_formula <- user_formula
                              self$inla_fitted_formula <- inla_fitted_formula
                              
                              self$inla_incl_spatial_term <- inla_incl_spatial_term
                              self$inla_stack <- inla_stack
                              self$inla_mesh <- inla_mesh
                              self$inla_mesh_settings <- inla_mesh_settings
                              
                              # Incl processing flags from hib_db
                              self$processing_flags <- db_obj$processing_flags
                              
                              # Calc n sites/species used in  model - need to calc spp with div metrics!!
                              self$processing_flags$N_spp <- n_spp
                              
                              if ("SSBS" %in% colnames(mod_df)){
                                self$processing_flags$N_sites <- dplyr::n_distinct(mod_df$SSBS)
                              } else if ("SSBS1" %in% colnames(mod_df)){
                                self$processing_flags$N_sites <- length(unique(c(mod_df$SSBS1, mod_df$SSBS2)))
                              }
                              
                              self$version <- db_obj$version
                              self$last_update <- lubridate::today()
                              
                              # %%%%%%%%%%%%%%%
                              
                            },
                            # Functions
                            
                            # Base printout function
                            print = function(){
                              stats <- list(number_sites = self$processing_flags$N_sites,
                                            number_species = self$processing_flags$N_spp ) 
                              # Metrics of fit?
                              m_summ <- summary(self$model)
                              
                              fixef_df <- m_summ$fixed %>%
                                as.data.frame(.) %>%
                                tibble::rownames_to_column("Coefficient") %>%
                                dplyr::tibble(.) %>%
                                dplyr::rename(q0.025 = `0.025quant`,
                                              q0.975 = `0.975quant`) %>%
                                dplyr::select(Coefficient, mean, sd, q0.025, q0.975) 
                              
                              ranef_chr <- purrr::map2(.x = m_summ$random.names,
                                                .y = m_summ$random.model,
                                                function(.x, .y){paste0("\t", .x, "\t ", .y)}) %>%
                                      unlist() %>%
                                      paste0(., collapse = "\n")
                              
                              
                              message(
                                paste0("Hibiiscus fitted INLA model\n\n",
                                       "(Last updated: ", crayon::green( self$last_update ), ")",
                                       "\n",
                                       "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                                       "\t Number of species:             ", crayon::green(stats$number_species), "\n",
                                       "\n\n",
                                       "--------------------\n\n",
                                       "Summary of fitted model:\n\n",
                                       "Fixed effects:"
                                       )
                                )
                              print(fixef_df)

                              message(paste0("\n\n",
                                             # "--------------------\n",
                                             "Random effects:\n",
                                             ranef_chr, "\n\n",
                                             "Performance metrics:\n",
                                             "\tDIC \t ", m_summ$dic$dic, "\n",
                                             "\twAIC \t ", m_summ$waic$waic, "\n"))
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
                            plot = function(type = "fixed-effects"){
                              assertthat::assert_that(is.character(type))
                              # stop("Function in progress...")
                              # Plot options
                              type = match.arg(type,
                                               c("fixed-effects",
                                                 # "fixed-effects posteriors",
                                                 "random",
                                                 "spatial-random-field"))
                                              
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
                                fixef_df <- self$model$summary.fixed
                                fixef_df$Coefficient <- rownames(fixef_df)
                                
                                # Set a baseline- Intercept
                                # %%%
                                int_id <- which(grepl("intercept", fixef_df$Coefficient, ignore.case = T))
                                # TODO: What if another covar also contains intercept...
                                if (length(int_id) > 1){
                                  message(crayon::yellow("Warning: Multiple 'Intercept' coefficients found, selecting the first...\n"))
                                  int_id <- int_id[1]
                                }
                                # If no explicit intercept, set as 1st value...
                                if (length(int_id) == 0){
                                  int_id <- 1
                                }
                                fixef_df$Coefficient <- factor(fixef_df$Coefficient)
                                fixef_df$Coefficient <- relevel(fixef_df$Coefficient,
                                                                ref = as.character(fixef_df$Coefficient[int_id]))
                                # %%%
                                plt <- ggplot(fixef_df) +
                                  geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                                  geom_errorbar(aes(x = Coefficient, 
                                                    ymin = `0.025quant`,
                                                    ymax = `0.975quant`),
                                                width = 0) +
                                  geom_point(aes(x = Coefficient,
                                                 y = mean)) +
                                  
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
                              if (type == "random"){
                                # if (is.null(self$user_formula$random)){
                                #   
                                # }
                                # Get random vars from random form
                                rand_vars <- all.vars(self$user_formula$random)
                                # Get random vars coefs 
                                rand_df <- self$model$summary.random[rand_vars] %>%
                                  dplyr::bind_rows(., .id = "Covar") %>%
                                  # order from small to large?
                                  dplyr::group_by(Covar) %>%
                                  dplyr::arrange(mean, .by_group = T) %>%
                                  dplyr::mutate(ID_ord = 1:n()) %>%
                                  dplyr::ungroup() %>%
                                  dplyr::mutate(sig_col = case_when(`0.025quant` > 0 ~ "+ve",
                                                                    `0.975quant` < 0 ~ "-ve",
                                                                    T ~ "n.s."))
                                
                                # colour as +ve or -ve 
                                # plot - no names...
                                plt <- ggplot(rand_df) +
                                  geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                                  geom_errorbar(aes(x = ID_ord, 
                                                    ymin = `0.025quant`,
                                                    ymax = `0.975quant`,
                                                    colour = sig_col),
                                                show.legend = F,
                                                alpha = 0.5,
                                                width = 0) +
                                  geom_point(aes(x = ID_ord,
                                                 y = mean,
                                                 colour = sig_col),
                                             show.legend = F,
                                             alpha = 0.7) +
                                  labs(y = "Estimate",
                                       x = element_blank(),
                                       title = "Random effect estimates") +
                                  scale_color_manual(values = c("+ve" = "dodgerblue3",
                                                                "-ve" = "firebrick3",
                                                                "n.s."= "grey30")) +
                                  facet_wrap(~Covar, scales = "free_x") +
                                  theme_bw() +
                                  basic_thm +
                                  theme(axis.text.x = element_blank(),
                                        aspect.ratio = 1)
                                
                                plot(plt)
                              }
                              # 3. spatial field
                              if (type == "spatial-random-field"){
                                # check spatial randome effect is True
                                if (!self$inla_incl_spatial_term){
                                  message(crayon::red("Error: plotting of spatial-random-field is only available for INLA models with a spatial random effect...\n"))     
                                  stop()
                                }
                                pix <- inlabru::pixels(self$inla_mesh)
                                proj <- inla.mesh.projector(self$inla_mesh, coordinates(pix))
                                pix$z <- inla.mesh.project(proj, self$model$summary.random$s$mean)
                                
                                plt <- ggplot() + 
                                  inlabru::gg(self$inla_mesh) +
                                  inlabru::gg(pix, aes(fill = z),
                                              alpha = 0.9) +
                                  scale_fill_gradient2(low = "firebrick3",
                                                       mid = "gainsboro",
                                                       high = "dodgerblue3",
                                                       midpoint = 0) + 
                                  labs(fill = "Spatial\nrandom\neffect") +
                                  coord_fixed() +
                                  theme_bw()
                                plot(plt)
                              }
                              # 4. priors
                              # 5. partial effects
                              
                              
                            }
                          ),
                          private = list()
)

