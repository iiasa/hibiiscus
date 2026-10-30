#' Template object for INLA model 
#' @title HibModelINLA
#' @description 
#' This is the template object for an [R6][R6::R6Class] `INLA` model class.
#' It contains functions for summarising and visualising the model results. 
#' @keywords classes
#' @import R6
#' @rdname HibModelINLA-Class
#' @export
HibModelINLA <- R6Class("HibModelINLA", 
                        
                          # Public data objects and functions
                          public = list(
                            #' @field dataset_name (`character`)\cr
                            #' Identifies the dataset used.
                            dataset_name = NULL,
                            
                            #' @field last_update (`character`)\cr
                            #' Identifies the most recent update to the object (YYYY-MM-DD).
                            last_update = NULL,
                            
                            #' @field data (`data.frame`)\cr
                            #' The data used in model fitting.
                            data = NULL,
                            
                            #' @field model (`inla`)\cr
                            #' The `INLA::inla` model object.
                            model = NULL,
                            
                            #' @field processing_flags (`list()`)\cr
                            #' Identifies the processing steps used for the HibDatabase 
                            #' object that was used for model fitting.
                            processing_flags = NULL,
                            
                            #' @field fixedef_coefficients (`data.frame`)\cr
                            #' The fixed effects coefficients.
                            fixedef_coefficients = NULL,
                            
                            #' @field user_formula (`formula` or `list`)\cr
                            #' The model formula supplied by the user (includes
                            #' the fixed effects and non-spatial random effects).
                            user_formula = NULL,
                            
                            #' @field inla_fitted_formula (`formula`)\cr
                            #' The fitted model formula (includes all effects).
                            inla_fitted_formula = NULL,
                            
                            #' @field inla_incl_spatial_term (`logical`)\cr
                            #' TRUE or FALSE, does the model include a spatial term.
                            inla_incl_spatial_term = NULL,
                            
                            #' @field inla_stack (`inla.data.stack`)\cr
                            #' The data.stack used in model fitting, an INLA object.
                            inla_stack = NULL,
                            
                            #' @field inla_mesh (`inla.mesh` or `NULL`)\cr
                            #' If used, the mesh used in model fitting.
                            inla_mesh = NULL,
                            
                            #' @field inla_mesh_settings (`list`)\cr
                            #' A list of specifications for generating the mesh,
                            #' if used in modelling. 
                            inla_mesh_settings = NULL,
                            
                            #' @field inla_mesh_source (`character` or `NULL`)\cr
                            #' A description of the mesh source - if not null, 
                            #' can be either "user_provided" or "mesh_settings".
                            inla_mesh_source = NULL,
                            
                            #' @field site_locations (`data.frame`)\cr
                            #' Contains columns for "SSBS", "Longitude" ,"Latitude".
                            site_locations = NULL,
                            
                            
                            # what else to include
                            # 
                            # inla output is fairly basic, but extensive...
                            # important to store random effects params/hypers
                            # need to store A matrix
                            # store priors etc
                            
                            
                            # Executed when created
                            #' @description
                            #' Creates a new instance of this [R6][R6::R6Class] class.
                            #' Note, this object is constructed via the [hib_fit_inla]
                            #' function found within this package.
                            #' @param db_obj (`HibDatabase`, or `list`)\cr
                            #'   The HibDatabase object/s used to fit the model.
                            #' @param mod_df (`data.frame`)\cr
                            #'   The data.frame of the data used in model fitting.
                            #' @param inla_obj (`inla`)\cr
                            #'   The inla model object.
                            #' @param fixef_df (`data.frame`)\cr
                            #'   The fixed effects from the fitted model.
                            #' @param metr_ls (`list`)\cr
                            #'   The list of model performance metrics.
                            #' @param user_formula (`formula` or `list`)\cr
                            #'   The model formula supplied by the user.
                            #' @param inla_fitted_formula (`formula`)\cr
                            #'   The fitted model formula.
                            #' @param inla_incl_spatial_term (`logical`)\cr
                            #'   TRUE or FALSE, does the model include a spatial term.
                            #' @param inla_stack (`inla.data.stack`)\cr
                            #'   The data.stack used in model fitting, an INLA object.
                            #' @param inla_mesh (`inla.mesh` or `NULL`)\cr
                            #'   If used, the mesh used in model fitting.
                            #' @param inla_mesh_settings (`list`)\cr
                            #'   A list of specifications for generating the mesh. 
                            #' @param inla_mesh_source (`character` or `NULL`)\cr
                            #'   A description of the mesh source.
                            #' @param n_spp (`integer`)\cr
                            #'   Number of species contributing to the fitted model.
                            #' @param n_taxa (`integer`)\cr
                            #'   Number of orignal taxa names contributing to the fitted model.
                            #' @param n_sites (`integer`)\cr
                            #'   Number of sites contributing to the fitted model.
                            #'   
                            #'   @return (`HibModelINLA`)
                            initialize = function(db_obj, 
                                                  mod_df,
                                                  inla_obj,
                                                  fixef_df,
                                                  user_formula,
                                                  inla_fitted_formula,
                                                  inla_incl_spatial_term,
                                                  inla_stack,
                                                  inla_mesh,
                                                  inla_mesh_settings,
                                                  inla_mesh_source,
                                                  n_spp,
                                                  n_taxa,
                                                  n_sites){

                              # Package check
                              if (!requireNamespace("INLA", quietly = TRUE)) {
                                stop(
                                  "Package \"INLA\" must be installed to use this class",
                                  call. = FALSE
                                )
                              }
                              # OR
                              # rlang::check_installed("INLA", reason = "to use INLA-based model fitting.")
                              
                              # - some assertion checks
                              assertthat::assert_that(
                                # Check database obj
                                inherits(db_obj, "HibDatabase"),
                                # Check model data
                                is.data.frame(mod_df),
                                # Check model
                                class(inla_obj) == "inla",
                                # Check coefficient df 
                                is.data.frame(fixef_df),
                                # Check formulas
                                class(user_formula) == "formula" || is.list(user_formula),
                                class(inla_fitted_formula) == "formula",
                                # Check inla specific data
                                inherits(inla_stack, "inla.data.stack"),
                                is.logical(inla_incl_spatial_term),
                                inherits(inla_mesh,'inla.mesh') || is.null(inla_mesh),
                                
                                is.list(inla_mesh_settings),
                                is.character(inla_mesh_source) || is.null(inla_mesh_source),
                                is.numeric(n_spp),
                                is.numeric(n_taxa),
                                is.numeric(n_sites)
                              )
                              
                              # ------ #
                              message(crayon::green("\nInitialising HibModelINLA object ...") )
                              
                              # %%%%%%%%%%%%%%%
                              
                              # Save inputs in self
                              # Set mod_df as data
                              self$data <- mod_df
                              # Set model as inla_obj
                              self$model <- inla_obj
                              # Set fixed_coefficients as fixef_df
                              self$fixedef_coefficients <- fixef_df
                              
                              self$user_formula <- user_formula
                              self$inla_fitted_formula <- inla_fitted_formula
                              
                              self$inla_incl_spatial_term <- inla_incl_spatial_term
                              self$inla_stack <- inla_stack
                              self$inla_mesh <- inla_mesh
                              self$inla_mesh_settings <- inla_mesh_settings
                              
                              self$inla_mesh_source <- inla_mesh_source
                              
                              # Incl processing flags from hib_db
                              self$processing_flags <- db_obj$processing_flags
                              self$processing_flags$N_spp <- n_spp
                              self$processing_flags$N_taxa <- n_taxa
                              self$processing_flags$N_sites <- n_sites
                            
                              self$dataset_name <- db_obj$dataset_name

                              # Store site locations
                              loc_df <- db_obj$sites[, c("SSBS", "Longitude" ,"Latitude")]
                              if ("SSBS" %in% colnames(mod_df)){
                                loc_df <- loc_df |> dplyr::filter(SSBS %in% mod_df$SSBS)
                              } else if ("SSBS1" %in% colnames(mod_df)){
                                loc_df <- loc_df |> dplyr::filter(SSBS %in% unique(c(mod_df$SSBS1, mod_df$SSBS2)))
                              }
                              self$site_locations <- loc_df
                              
                              self$last_update <- lubridate::now()
                              
                              # %%%%%%%%%%%%%%%
                              
                            },
                            # Functions
                            
                            #' @description
                            #' Base printout function presents core information about the 
                            #' object.
                            #'
                            #' @return (`NULL`).
                            print = function(){
                              stats <- list(number_sites = self$processing_flags$N_sites,
                                            number_taxa = self$processing_flags$N_taxa,
                                            number_species = self$processing_flags$N_spp ) 
                              # Metrics of fit
                              m_summ <- summary(self$model)
                              
                              # neatly print formula
                              form_chr <- as.character(self$user_formula$full)
                              form_chr <- paste(form_chr[2], form_chr[1], form_chr[3])
                              if (self$inla_incl_spatial_term){
                                form_chr <- paste(form_chr,
                                                  "+ f(s, model = site.spde)")
                              }
                              
                              fixef_df <- self$fixedef_coefficients |>
                                dplyr::select(Coefficient_name, Coefficient_value,
                                              sd, loCI, hiCI)
                              
                              ranef_chr <- purrr::map2(.x = m_summ$random.names,
                                                .y = m_summ$random.model,
                                                function(.x, .y){paste0("\t", .x, "\t ", .y)}) |>
                                      unlist() |>
                                      paste0(collapse = "\n")
                              
                              
                              message(
                                paste0("Hibiiscus fitted model (INLA)\n\n",
                                       "Last updated: ", crayon::green( self$last_update ), 
                                       "\n",
                                       # TODO: Dataset name/s ?
                                       "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                                       "\t Number of species:             ", crayon::green(stats$number_species), "\n",
                                       "\t Number of taxa (orig. names):  ", crayon::green(stats$number_taxa), "\n",
                                       "\n",
                                       "--------------------\n\n",
                                       "Summary of fitted model\n\n",
                                       "Formula:\n",
                                       crayon::green(form_chr), "\n\n",
                                       "Fixed effects estimates:\n"
                                       )
                                )
                              print(fixef_df)

                              message(paste0("\n\n",
                                             # "--------------------\n",
                                             # "Random effects:\n",
                                             # ranef_chr, "\n\n",
                                             "Performance metrics:\n",
                                             "\tDIC \t ", m_summ$dic$dic, "\n",
                                             "\twAIC \t ", m_summ$waic$waic, "\n"))
                            },
                            
                            
                            #' @description
                            #' Function to generate figures summarising the model.
                            #' 
                            #' @param type (`character()`) Options are: "fixed-effects",
                            #' "random-effects", "spatial-random-field", "map-sites".
                            #' @return (`ggplot`).
                            plot = function(type = "fixed-effects"){
                              assertthat::assert_that(is.character(type))
                              # Plot options
                              # TODO: covar-distributions
                              type = match.arg(type,
                                               c("fixed-effects",
                                                 # "fixed-effects posteriors",
                                                 "random-effects",
                                                 "spatial-random-field",
                                                 "map-sites"))
                                              
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
                              
                              # 1. fixed effects - summary
                              if (type == "fixed-effects"){
                                # use fixedef_coefficients df
                                fixef_df <- self$fixedef_coefficients
                                
                                # Convert chr col of Coefficient_name to factor, based on order of terms in df
                                fixef_df$Coefficient_name <- factor(fixef_df$Coefficient_name,
                                                                    levels = fixef_df$Coefficient_name)


                                # Use Coefficient_name as x
                                plt <- ggplot2::ggplot(fixef_df) +
                                  ggplot2::geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                                  ggplot2::geom_errorbar(ggplot2::aes(x = Coefficient_name,
                                                ymin = loCI,
                                                ymax = hiCI),
                                            width = 0) +
                                  ggplot2::geom_point(ggplot2::aes(x = Coefficient_name,
                                                                   y = Coefficient_value)) +
                                  
                                  ggplot2::labs(y = "Estimate",
                                       x = ggplot2::element_blank(),
                                       title = "Fixed effect estimates") +
                                  ggplot2::theme_bw() +
                                  basic_thm +
                                  ggplot2::coord_flip()
                                
                              }
                              #    fixed effects - posteriors
                              #    TODO
                              # 2. random iid
                              if (type == "random-effects"){
                                
                                # Get random vars from random form
                                rand_vars <- all.vars(self$user_formula$random)
                                # Get random vars coefs 
                                
                                # Check inla summary.random - what if random slopes??
                                rand_df <- self$model$summary.random[rand_vars] |>
                                  dplyr::bind_rows(.id = "Covar") |>
                                  # order from small to large
                                  dplyr::group_by(Covar) |>
                                  dplyr::arrange(mean, .by_group = T) |>
                                  dplyr::mutate(ID_ord = 1:dplyr::n()) |>
                                  dplyr::ungroup() |>
                                  dplyr::mutate(sig_col = dplyr::case_when(`0.025quant` > 0 ~ "+ve",
                                                                    `0.975quant` < 0 ~ "-ve",
                                                                    T ~ "n.s."))
                                
                                # colour as +ve or -ve 
                                # plot - no names...
                                plt <- ggplot2::ggplot(rand_df) +
                                  ggplot2::geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                                  ggplot2::geom_errorbar(ggplot2::aes(x = ID_ord, 
                                                    ymin = `0.025quant`,
                                                    ymax = `0.975quant`,
                                                    colour = sig_col),
                                                show.legend = F,
                                                alpha = 0.5,
                                                width = 0) +
                                  ggplot2::geom_point(ggplot2::aes(x = ID_ord,
                                                 y = mean,
                                                 colour = sig_col),
                                             show.legend = F,
                                             alpha = 0.7) +
                                  ggplot2::labs(y = "Estimate",
                                       x = ggplot2::element_blank(),
                                       title = "Random effect estimates") +
                                  ggplot2::scale_color_manual(values = c("+ve" = "dodgerblue3",
                                                                "-ve" = "firebrick3",
                                                                "n.s."= "grey30")) +
                                  ggplot2::facet_wrap(~Covar, scales = "free_x") +
                                  ggplot2::theme_bw() +
                                  basic_thm +
                                  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                        aspect.ratio = 1)
                                
                              }
                              # 3. spatial field
                              if (type == "spatial-random-field"){
                                # check spatial randome effect is True
                                if (!self$inla_incl_spatial_term){
                                  message(crayon::red("Error: plotting of spatial-random-field is only available for INLA models with a spatial random effect...\n"))     
                                  stop()
                                }
                                
                                # Facilitate for 3d meshes too
                                # %%%%%
                                if (self$inla_mesh_settings$mesh_dim == "3d"){
                                  gproj <- INLA::inla.mesh.projector(self$inla_mesh,  dims = c(360, 180))
                                  smean <- INLA::inla.mesh.project(gproj,
                                                                   self$model$summary.random$s$mean)
                                  
                                  # TODO: better plotting...
                                  
                                  # Convert smean to dataframe
                                  coord <- expand.grid(x = 1:360, y = 1:180) 
                                  out <- transform(coord, z = smean[as.matrix(coord)])
                                  out$x <- out$x - 180.5
                                  out$y <- out$y - 90.5
                                  
                                  ggplot2::ggplot(out) +
                                    ggplot2::geom_tile(ggplot2::aes(x = x,
                                                                    y = y,
                                                                    fill = z),
                                                       alpha = 0.9) +
                                    # TODO: add locations of data?
                                    # TODO: add underlying geogr map of country outlines
                                    # ggplot2::geom_point(ggplot2::aes(x = Longitude, 
                                    #                                  y = Latitude),
                                    #                     alpha = 0.5,
                                    #                     data = ta_inla_s$data) +
                                    ggplot2::scale_fill_gradient2(low = "firebrick3",
                                                                  mid = "gainsboro",
                                                                  high = "dodgerblue3",
                                                                  midpoint = 0) +
                                    ggplot2::labs(fill = "Spatial\nrandom\neffect") +
                                    ggplot2::coord_fixed(expand = F) +
                                    ggplot2::theme_bw() +
                                    basic_thm
                                  
                                }
                                
                                
                                # %%%%%
                                # If inla mesh is 2d...
                                if (self$inla_mesh_settings$mesh_dim == "2d"){
                                  # TODO: fix dependencies here!
                                  if (!requireNamespace("inlabru", quietly = TRUE)) {
                                    stop(
                                      "Package \"inlabru\" must be installed to plot 2d random fields...",
                                      call. = FALSE
                                    )
                                  }
                                  if (!requireNamespace("sp", quietly = TRUE)) {
                                    stop(
                                      "Package \"sp\" must be installed to plot 2d random fields...",
                                      call. = FALSE
                                    )
                                  }
                                  
                                  pix <- inlabru::pixels(self$inla_mesh)
                                  # TODO: conv sp:: to sf::st_coordinates... - need to convert pix first... st_as_sf() ?
                                  proj <- INLA::inla.mesh.projector(self$inla_mesh, sp::coordinates(pix))
                                  pix$z <- INLA::inla.mesh.project(proj, self$model$summary.random$s$mean)
                                  
                                  plt <- ggplot2::ggplot() + 
                                    inlabru::gg(self$inla_mesh) +
                                    inlabru::gg(pix, ggplot2::aes(fill = z),
                                                alpha = 0.9) +
                                    ggplot2::scale_fill_gradient2(low = "firebrick3",
                                                                  mid = "gainsboro",
                                                                  high = "dodgerblue3",
                                                                  midpoint = 0) + 
                                    ggplot2::labs(fill = "Spatial\nrandom\neffect") +
                                    ggplot2::coord_fixed() +
                                    ggplot2::theme_bw()
                                }
                                
                                
                              }
                              # 4. priors
                              # 5. partial effects
                              
                              if (type == "map-sites"){
                                
                                # Just get a coarse-ish world map...
                                tmp <- rnaturalearth::ne_countries(scale = 50, returnclass = 'sf')
                                
                                resp_var <- all.vars(self$user_formula$response_fixed[[2]])
                                resp_lab <- gsub("_", "\n", resp_var)
                                
                                loc_df <- self$site_locations
                                # %%%
                                
                                plt <- ggplot2::ggplot() +
                                  ggplot2::geom_sf(data = tmp,
                                                   fill = "lightgrey",
                                                   alpha = 0.4,
                                                   show.legend = F) +
                                  ggplot2::theme_bw() +
                                  basic_thm 
                                # %%%
                                
                                plt_df <- self$data
                                
                                if ("SSBS" %in% colnames(plt_df)){
                                  plt_df <- plt_df |>
                                    dplyr::left_join(loc_df, by = "SSBS")
                                  
                                }
                                if ("SSBS1" %in% colnames(plt_df)){
                                  plt_df <- plt_df |>
                                    # get the lat and long for s1
                                    dplyr::left_join(loc_df, by = c("SSBS1" = "SSBS")) |>
                                    dplyr::rename(SSBS1_latitude = Latitude,
                                                  SSBS1_longitude = Longitude) |>
                                    # get the lat and long for s2
                                    dplyr::left_join(loc_df, by = c("SSBS2" = "SSBS")) |>
                                    dplyr::rename(SSBS2_latitude = Latitude,
                                                  SSBS2_longitude = Longitude) |>
                                    # calculate the geographic midpoint between s1 and s2 
                                    dplyr::mutate(Latitude = (SSBS1_latitude + SSBS2_latitude) / 2,
                                                  Longitude = (SSBS1_longitude + SSBS2_longitude) / 2)
                                  
                                  # Or plot line segment between s1 and s2...
                                }
                                
                                # Plot site locations - colour as response variable - untrans...
                                
                                # Get x/y limits for sites...
                                x_rng <- range(plt_df$Longitude, na.rm = T)
                                y_rng <- range(plt_df$Latitude, na.rm = T)
                                
                                
                                # buffer by 5 degrees
                                x_rng[1] <- round(x_rng[1] - 5, 2)
                                x_rng[2] <- round(x_rng[2] + 5, 2)
                                
                                y_rng[1] <- round(y_rng[1] - 5, 2)
                                y_rng[2] <- round(y_rng[2] + 5, 2)
                                
                                x_rng[1] <- max(c(-180, x_rng[1]))
                                x_rng[2] <- min(c(180, x_rng[2]))
                                
                                y_rng[1] <- max(c(-90, y_rng[1]))
                                y_rng[2] <- min(c(90, y_rng[2]))
                                
                                
                                plt <- plt +
                                  ggplot2::geom_point(ggplot2::aes(x = Longitude,
                                                                   y = Latitude,
                                                                   colour = !!dplyr::sym(resp_var)),
                                                      alpha = 0.7,
                                                      data = plt_df) +
                                  ggplot2::labs(colour = resp_lab,
                                                x = ggplot2::element_blank(),
                                                y = ggplot2::element_blank(),
                                                title = "Geographic distribution of sites in model") +
                                  ggplot2::scale_color_gradient(low = "firebrick3",
                                                                high = "dodgerblue3") +
                                  ggplot2::xlim(x_rng) +
                                  ggplot2::ylim(y_rng) 
                                
                              }
                              
                              (plt)
                            }
                            
                            # TODO: check model performance
                          ),
                          private = list()
)

