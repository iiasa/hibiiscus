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
                      ci_type = NULL,
                      metrics = NULL,
                      processing_flags = NULL,
                      site_locations = NULL,
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
                                            ci_type,
                                            metr_vec,
                                            n_spp
                                            # set_up_specs
                                            # spp df?
                                            ){
                        # - some assertion checks
                        assertthat::assert_that(
                          # Check database obj
                          # db_obj$checkself(),
                          inherits(db_obj, "Database") || is.list(db_obj),
                          # Check model data
                          is.data.frame(mod_df),
                          # Check model
                          assertthat::assert_that(class(lme4_obj) == "lmerMod" || class(lme4_obj) == "glmerMod"), 
                          # Check coefficient df and metrics
                          is.data.frame(fixef_df),
                          is.character(ci_type) || is.na(ci_type),
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
                        
                        self$ci_type <- ci_type
                        
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
                        
                        if (inherits(db_obj, "Database")){
                          self$version <- db_obj$version  
                        }
                        if (is.list(db_obj)){
                          self$version <- lapply(db_obj, function(x){x$version})
                        }
                        

                        # # TODO: store site locations
                        loc_df <- db_obj$sites[, c("SSBS", "Longitude" ,"Latitude")]
                        if ("SSBS" %in% colnames(mod_df)){
                          loc_df <- loc_df %>% dplyr::filter(SSBS %in% mod_df$SSBS)
                        } else if ("SSBS1" %in% colnames(mod_df)){
                          loc_df <- loc_df %>% dplyr::filter(SSBS %in% unique(c(mod_df$SSBS1, mod_df$SSBS2)))
                        }
                        self$site_locations <- loc_df

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
                        # N spp/sites...
                        # print(summary(self$model))
                        stats <- list(number_sites = self$processing_flags$N_sites,
                                      # number_lu = length( unique(self$data$Predominant_land_use) ),
                                      # number_lui = length(unique(self$data$Use_intensity)),
                                      number_species = paste(self$processing_flags$N_spp, collapse = ", ") ) #,
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
                                 "--------------------\n\n",
                                 "Summary of fitted model\n"
                                 # "Fixed effects:"
                          )
                        )
                        # print(self$fixedef_coefficients)
                        print(summary(self$model))
                        metr <- pad_fn(self$metrics, 3, F)
                        names(metr) <- names(self$metrics)
                        
                        message(paste0("\n\n",
                                       # "--------------------\n",
                                       # "Random effects:\n",
                                       # ranef_chr, "\n\n",
                                       "Performance metrics:\n",
                                       "\tR-squared (marg.) \t ", metr["R2m"], "\n",
                                       "\tR-squared (cond.) \t ", metr["R2c"], "\n",
                                       "\tAIC               \t ", metr["AIC"], "\n",
                                       "\tAICc              \t ", metr["AICc"], "\n",
                                       "\tBIC               \t ", metr["BIC"], "\n"))
                        
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
                          fixef_df$Coefficient_name <- factor(fixef_df$Coefficient_name)
                          if ("(Intercept)" %in% fixef_df$Coefficient_name){
                            fixef_df$Coefficient <- relevel(fixef_df$Coefficient_name,
                                                            ref = "(Intercept)")
                          }
                          
                          # %%%
                          plt <- ggplot2::ggplot(fixef_df) +
                            ggplot2::geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                            ggplot2::geom_point(ggplot2::aes(x = Coefficient_name,
                                           y = Coefficient_value)) +
                            
                            ggplot2::labs(y = "Estimate",
                                 x = ggplot2::element_blank(),
                                 title = "Fixed effect estimates") +
                            ggplot2::theme_bw() +
                            basic_thm +
                            ggplot2::coord_flip()
                          
                          if (!is.na(self$ci_type)){
                            plt <- plt +
                              ggplot2::geom_errorbar(ggplot2::aes(x = Coefficient_name,
                                                ymin = loCI,
                                                ymax = hiCI),
                                            width = 0) 
                          }
                          
                        }
                        #    fixed effects - posteriors
                        # 2. random iid
                        if (type == "random-effects"){
                         
                          rand_vars <- all.vars(self$user_formula$random)
                          
                          rand_df <- as.data.frame(lme4::ranef(self$model)) %>%
                            dplyr::mutate(lwr = condval - 1.96*condsd, 
                                          upr = condval + 1.96*condsd) %>%
                            # order from small to large?
                            dplyr::group_by(grpvar, term) %>%
                            dplyr::arrange(condval, .by_group = T) %>%
                            dplyr::mutate(ID_ord = 1:dplyr::n()) %>%
                            dplyr::ungroup() %>%
                            dplyr::mutate(sig_col = dplyr::case_when(lwr > 0 ~ "+ve",
                                                              upr < 0 ~ "-ve",
                                                              T ~ "n.s."))
                          
                          # colour as +ve or -ve 
                          # plot - no names...
                          plt <- ggplot2::ggplot(rand_df) +
                            ggplot2::geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                            ggplot2::geom_errorbar(ggplot2::aes(x = ID_ord, 
                                              ymin = lwr,
                                              ymax = upr,
                                              colour = sig_col),
                                          show.legend = F,
                                          alpha = 0.5,
                                          width = 0) +
                            ggplot2::geom_point(ggplot2::aes(x = ID_ord,
                                           y = condval,
                                           colour = sig_col),
                                       show.legend = F,
                                       alpha = 0.7) +
                            ggplot2::labs(y = "Estimate",
                                 x = ggplot2::element_blank(),
                                 title = "Random effect estimates") +
                            ggplot2::scale_color_manual(values = c("+ve" = "dodgerblue3",
                                                          "-ve" = "firebrick3",
                                                          "n.s."= "grey30")) +
                            ggplot2::facet_grid(term~grpvar, scales = "free_x") +
                            ggplot2::theme_bw() +
                            basic_thm +
                            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                  aspect.ratio = 1)
                          
                        }
                        
                        if (type == "map-sites"){
                          
                          
                          # Just get a coarse-ish world map...
                          tmp <- rnaturalearth::ne_countries(scale = 50, returnclass = 'sf')
                          # NOTE: This differe per model type!!
                          resp_var <- all.vars(self$model@call$formula[[2]])
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
                            plt_df <- plt_df %>%
                              dplyr::left_join(loc_df, by = "SSBS")
                            
                          }
                          if ("SSBS1" %in% colnames(plt_df)){
                            plt_df <- plt_df %>%
                              # get the lat and long for s1
                              dplyr::left_join(loc_df, by = c("SSBS1" = "SSBS")) %>%
                              dplyr::rename(SSBS1_latitude = Latitude,
                                            SSBS1_longitude = Longitude) %>%
                              # get the lat and long for s2
                              dplyr::left_join(loc_df, by = c("SSBS2" = "SSBS")) %>%
                              dplyr::rename(SSBS2_latitude = Latitude,
                                            SSBS2_longitude = Longitude) %>%
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
                                          y = ggplot2::element_blank()) +
                            ggplot2::scale_color_gradient(low = "firebrick3",
                                                          high = "dodgerblue3") +
                            ggplot2::xlim(x_rng) +
                            ggplot2::ylim(y_rng) 
                          
                        }
                        
                        (plt)
                      }
                    ),
                    private = list()
)