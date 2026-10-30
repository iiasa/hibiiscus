#' Template object for brms model 
#' @title HibModelBRMS
#' @description 
#' This is the template object for an [R6][R6::R6Class] `brms` model class.
#' It contains functions for summarising and visualising the model results. 
#' @keywords classes
#' @import R6
#' @rdname HibModelBRMS-Class
#' @export
HibModelBRMS <- R6Class("HibModelBRMS",
                          
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
                            
                            #' @field model (`brmsfit`)\cr
                            #' The brms model object.
                            model = NULL,
                            
                            #' @field extra_args (`list()`)\cr
                            #' Identifies addition arguments supplied during model 
                            #' fitting (found in '...' when calling `hib_fit_brms`).
                            extra_args = NULL,
                            
                            #' @field fixedef_coefficients (`data.frame`)\cr
                            #' The fixed effects coefficients.
                            fixedef_coefficients = NULL,
                            
                            #' @field metrics (`list`)\cr
                            #' List containing model performance metrics.
                            metrics = NULL,
                            
                            #' @field processing_flags (`list()`)\cr
                            #' Identifies the processing steps used for the HibDatabase 
                            #' object that was used for model fitting.
                            processing_flags = NULL,
                            
                            #' @field site_locations (`data.frame`)\cr
                            #' Contains columns for "SSBS", "Longitude" ,"Latitude"
                            site_locations = NULL,
                            
                            # TODO: what else to include
                            
                            # Executed when created
                            #' @description
                            #' Creates a new instance of this [R6][R6::R6Class] class.
                            #' Note, this object is constructed via the [hib_fit_brms]
                            #' function found within this package.
                            #' @param db_obj (`HibDatabase`, or `list`)\cr
                            #'   The HibDatabase object/s used to fit the model.
                            #' @param mod_df (`data.frame`)\cr
                            #'   The data.frame of the data used in model fitting.
                            #' @param brms_obj (`brmsfit`)\cr
                            #'   The brms model object.
                            #' @param fixef_df (`data.frame`)\cr
                            #'   The fixed effects from the fitted model.
                            #' @param metr_ls (`list`)\cr
                            #'   The list of model performance metrics.
                            #' @param n_spp (`integer` or `list`)\cr
                            #'   Number of species contributing to the fitted model.
                            #' @param n_taxa (`integer` or `list`)\cr
                            #'   Number of orignal taxa names contributing to the fitted model.
                            #' @param n_sites (`integer` or `list`)\cr
                            #'   Number of sites contributing to the fitted model.
                            #' @param extra_args (`list`)\cr
                            #'   Additional arguments specified during model fitting.
                            #'   
                            #'   @return (`HibModelBRMS`)
                            initialize = function(db_obj, 
                                                  mod_df,
                                                  brms_obj,
                                                  fixef_df,
                                                  metr_ls,
                                                  n_spp,
                                                  n_taxa,
                                                  n_sites,
                                                  extra_args){
                              # - some assertion checks
                              assertthat::assert_that(
                                # Check database obj
                                inherits(db_obj, "HibDatabase") ||  ( is.list(db_obj) & all(unlist(lapply(db_obj, function(x){inherits(x, "HibDatabase")}))) ), 
                                # Check model data
                                is.data.frame(mod_df),
                                # Check model
                                class(brms_obj) == "brmsfit",
                                # Check coefficient df and metrics
                                is.data.frame(fixef_df),
                                is.list(metr_ls),
                                is.numeric(n_spp) || ( is.list(n_spp) & all(unlist(lapply(n_spp, function(x){is.numeric(x)}))) ),
                                is.numeric(n_taxa) || ( is.list(n_taxa) & all(unlist(lapply(n_taxa, function(x){is.numeric(x)}))) ),
                                is.numeric(n_sites) || ( is.list(n_sites) & all(unlist(lapply(n_sites, function(x){is.numeric(x)}))) ),
                                is.list(extra_args)
                              )
                              
                              # ------ #
                              message(crayon::green("\nInitialising HibModelBRMS object ...") )
                              
                              # %%%%%%%%%%%%%%%
                              
                              # Save inputs in self
                              # Set mod_df as data
                              self$data <- mod_df
                              # Set model as brms_obj
                              self$model <- brms_obj

                              self$extra_args <- extra_args
                              
                              # Set fixed_coefficients as fixef_df
                              self$fixedef_coefficients <- fixef_df
                              # Set metrics as metr_ls
                              self$metrics <- metr_ls
                              
                              # Incl processing flags from hib_db
                              if (inherits(db_obj, "HibDatabase")){
                                self$processing_flags <- db_obj$processing_flags
                                self$processing_flags$N_spp <- n_spp
                                self$processing_flags$N_taxa <- n_taxa
                                self$processing_flags$N_sites <- n_sites
                              }
                              if (is.list(db_obj)){
                                self$processing_flags <- lapply(db_obj, function(x){x$processing_flags})
                                for (i in 1:length(db_obj)){
                                  self$processing_flags[[i]]$N_spp <- n_spp[[i]]
                                  self$processing_flags[[i]]$N_taxa <- n_taxa[[i]]
                                  self$processing_flags[[i]]$N_sites <- n_sites[[i]]
                                }
                              }
                            
                              # TODO: what about additional data (data2), e.g. covariance structures
                              
                              if (inherits(db_obj, "HibDatabase")){
                                self$dataset_name <- db_obj$dataset_name
                                loc_df <- db_obj$sites[, c("SSBS", "Longitude" ,"Latitude")]
                              }
                              if (is.list(db_obj)){
                                self$dataset_name <- lapply(db_obj, function(x){x$dataset_name})
                                loc_df <- lapply(db_obj, function(x){x$sites[, c("SSBS", "Longitude" ,"Latitude")]}) |>
                                  dplyr::bind_rows()
                              }
                            

                              # Store site locations
                              if ("SSBS" %in% colnames(mod_df)){
                                loc_df <- loc_df |> dplyr::filter(SSBS %in% mod_df$SSBS)
                              } else if ("SSBS1" %in% colnames(mod_df)){
                                loc_df <- loc_df |> dplyr::filter(SSBS %in% unique(c(mod_df$SSBS1, mod_df$SSBS2)))
                              }
                              self$site_locations <- loc_df

                              self$last_update <- Sys.time()
                              
                              # %%%%%%%%%%%%%%%
                           
                            },
                            # Functions
                            
                            #' @description
                            #' Base printout function presents core information about the 
                            #' object.
                            #'
                            #' @return (`NULL`).
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
                              # N spp/sites...
                              if (is.null(names(self$processing_flags))){
                                n_sites <- lapply(self$processing_flags, function(x){x$N_sites}) |>
                                  unlist() |>
                                  paste(collapse = ", ")
                                n_spp <- lapply(self$processing_flags, function(x){x$N_spp}) |>
                                  unlist() |>
                                  paste(collapse = ", ")
                                n_taxa <- lapply(self$processing_flags, function(x){x$N_taxa}) |>
                                  unlist() |>
                                  paste(collapse = ", ")
                                stats <- list(number_sites = n_sites,
                                              number_taxa = n_taxa,
                                              number_species = n_spp)
                              } else {
                                stats <- list(number_sites = self$processing_flags$N_sites,
                                              number_species = self$processing_flags$N_spp,
                                              number_taxa = self$processing_flags$N_taxa)
                              }

                              message(
                                paste0("Hibiiscus fitted model (brms)\n",
                                       "Last updated: ", crayon::green( self$last_update ),
                                       "\n",
                                       # TODO: Dataset name/s ?
                                       "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                                       # "\t Land use (class x intensity):  ", crayon::green(stats$number_lu), " x ",
                                       # crayon::green(stats$number_lui), "\n",
                                       "\t Number of species:             ", crayon::green(stats$number_species), "\n",
                                       "\t Number of taxa (orig. names):  ", crayon::green(stats$number_taxa), "\n",
                                       "\n",
                                       "--------------------\n\n",
                                       "Summary of fitted model\n"
                                       # "Fixed effects:"
                                )
                              )
                              print(summary(self$model))
                              
                              # Get loo/waic...
                              metr <- c(self$metrics$loo$estimates["looic","Estimate"],
                                        self$metrics$waic$estimates["waic","Estimate"])
                              metr <- pad_fn(metr, 3, F)
                              names(metr) <- c("looic", "waic")
                              
                              message(paste0("\n\n",
                                             # "--------------------\n",
                                             # "Random effects:\n",
                                             # ranef_chr, "\n\n",
                                             "Performance metrics:\n",
                                             "\tLOO IC \t ", metr["looic"], "\n",
                                             "\tWAIC   \t ", metr["waic"], "\n"))
                              
                            },
                          
                            
                            #' @description
                            #' Function to generate figures summarising the model.
                            #' 
                            #' @param type (`character()`) Options are: "fixed-effects",
                            #' "random-effects", "map-sites".
                            #' @return (`ggplot`).
                            plot = function(type = "fixed-effects"){
                              rlang::check_installed("ggplot2", reason = "to plot hibiiscus objects.")
                              
                              assertthat::assert_that(is.character(type))
                
                              # Plot options
                              type = match.arg(type,
                                               c("fixed-effects",
                                                 "random-effects",
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
                                fixef_df <- self$fixedef_coefficients
                                
                                # Set a baseline- Intercept
                                # Convert chr col of Coefficient_name to factor, based on order of terms in df
                                fixef_df$Coefficient_name <- factor(fixef_df$Coefficient_name,
                                                                    levels = fixef_df$Coefficient_name)
  
                                # %%%
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
                              # 2. random iid
                              if (type == "random-effects"){
                                rlang::check_installed(c("brms", "reshape2"),
                                                       reason = "to plot brms random effects.")
                                # TODO: ranef(m) produces a list of 3D arrays
                                # dim 1 = factor levels
                                # dim 2 = summary stats
                                # dim 3 = group-level effects
                                
                                # if random slopes, then x is 3d, not 2d
                                # coef nm, (intercept or LU) also in colnames
                                rand_df <- purrr::map(brms::ranef(self$model), 
                                           function(x){
                                             
                                             # melt - then colnames are coef_level, coef_stat, coef_term
                                             tmp <- reshape2::melt(x) |>
                                               # widen across stat
                                               tidyr::pivot_wider(names_from = "Var2",
                                                                  values_from = "value") |>
                                               dplyr::rename(Coef_level = Var1,
                                                             Coef_term = Var3,
                                                             Coefficient_value = Estimate,
                                                             Coefficient_error = Est.Error,
                                                             loCI = Q2.5,
                                                             hiCI = Q97.5)
                                            
                                             tmp
                                             }) |>
                                  dplyr::bind_rows(.id = "Coefficient_group") |>
                                  dplyr::group_by(Coefficient_group, Coef_term) |>
                                  dplyr::arrange(Coefficient_value, .by_group = T) |>
                                  dplyr::mutate(ID_ord = 1:dplyr::n()) |>
                                  dplyr::ungroup() |>
                                  dplyr::mutate(sig_col = dplyr::case_when(loCI > 0 ~ "+ve",
                                                                    hiCI < 0 ~ "-ve",
                                                                    T ~ "n.s."))
                                
                                
                                # colour as +ve or -ve 
                                # plot - no names...
                                plt <- ggplot2::ggplot(rand_df) +
                                  ggplot2::geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                                  ggplot2::geom_errorbar(ggplot2::aes(x = ID_ord, 
                                                    ymin = loCI,
                                                    ymax = hiCI,
                                                    colour = sig_col),
                                                show.legend = F,
                                                alpha = 0.5,
                                                width = 0) +
                                  ggplot2::geom_point(ggplot2::aes(x = ID_ord,
                                                 y = Coefficient_value,
                                                 colour = sig_col),
                                             show.legend = F,
                                             alpha = 0.7) +
                                  ggplot2::labs(y = "Estimate",
                                       x = ggplot2::element_blank(),
                                       title = "Random effect estimates") +
                                  ggplot2::scale_color_manual(values = c("+ve" = "dodgerblue3",
                                                                "-ve" = "firebrick3",
                                                                "n.s."= "grey30")) +
                                  ggplot2::facet_grid(Coef_term~Coefficient_group, scales = "free_x") +
                                  ggplot2::theme_bw() +
                                  basic_thm +
                                  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                        aspect.ratio = 1)

                              }
                              
                              if (type == "map-sites"){
                                
                                # Just get a coarse-ish world map...
                                rlang::check_installed(c("rnaturalearth", "rnaturalearthdata", "sf"),
                                                       reason = "to plot model locations.")
                                tmp <- rnaturalearth::ne_countries(scale = 50, returnclass = 'sf')
                                
                                resp_var <- all.vars(stats::formula(self$model)$formula[[2]])
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
                            },
                            
                            #' @description
                            #' Function to assess the model using the performance package.
                            #' 
                            #' @return (`NULL`).
                            check_model = function(){
                              message("Checking model...\n\tNote this function is still in development...")
                              m <- self$model
                              performance::check_model(m)
                            }
                          ),
                          private = list()
)