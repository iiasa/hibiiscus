#' Template object for lme4 model 
#' @title HibModelLME4
#' @description 
#' This is the template object for an [R6][R6::R6Class] lme4 model class.
#' It contains functions for summarising and visualising the model results. 
#' @keywords classes
#' @import R6
#' @rdname HibModelLME4-Class
#' @export
HibModelLME4 <- R6Class("HibModelLME4",# "hib_lme4_model",
                    
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
                      
                      #' @field model (`merMod`)\cr
                      #' The lme4 model object.
                      model = NULL,
                      
                      #' @field fixedef_coefficients (`data.frame`)\cr
                      #' The fixed effects coefficients.
                      fixedef_coefficients = NULL,
                      
                      #' @field ci_type (`character`, or `NA`)\cr
                      #' Identifies the types of confidence intervals calculated.
                      ci_type = NULL,
                      
                      #' @field metrics (`numeric`)\cr
                      #' Named vector containing model performance metrics.
                      metrics = NULL,
                      
                      #' @field processing_flags (`list()`)\cr
                      #' Identifies the processing steps used for the HibDatabase 
                      #' object that was used for model fitting.
                      processing_flags = NULL,
                      
                      #' @field site_locations (`data.frame`)\cr
                      #' Contains columns for "SSBS", "Longitude" ,"Latitude".
                      site_locations = NULL,
                      
                      # TODO: what else to include
                      # 
                      # lme4 output is pretty comprehensive, this is just a wrapper...
                      # need optimal ranef struct, if found
                      # need model ranks, if dredged..
                      
                      # Executed when created
                      #' @description
                      #' Creates a new instance of this [R6][R6::R6Class] class.
                      #' Note, this object is constructed via the [hib_fit_lme4]
                      #' function found within this package.
                      #' @param db_obj (`HibDatabase`, or `list`)\cr
                      #'   The HibDatabase object/s used to fit the model.
                      #' @param mod_df (`data.frame`)\cr
                      #'   The data.frame of the data used in model fitting.
                      #' @param lme4_obj (`merMod`)\cr
                      #'   The lme4 model object.
                      #' @param fixef_df (`data.frame`)\cr
                      #'   The fixed effects from the fitted model.
                      #' @param ci_type (`character`, or `NA`)\cr
                      #'   The type of confidence intervals calculated.
                      #' @param metr_vec (`numeric`)\cr
                      #'   The named vector of model performance metrics.
                      #' @param n_spp (`integer` or `list`)\cr
                      #'   Number of species contributing to the fitted model.
                      #' @param n_taxa (`integer` or `list`)\cr
                      #'   Number of orignal taxa names contributing to the fitted model.
                      #' @param n_sites (`integer` or `list`)\cr
                      #'   Number of sites contributing to the fitted model.
                      initialize = function(db_obj, # ??
                                            mod_df,
                                            lme4_obj,
                                            fixef_df,
                                            ci_type,
                                            metr_vec,
                                            n_spp,
                                            n_taxa,
                                            n_sites
                                            # set_up_specs
                                            # spp df?
                                            ){
                        # - some assertion checks
                        assertthat::assert_that(
                          # Check database obj
                          # db_obj$checkself(),
                          inherits(db_obj, "HibDatabase") || ( is.list(db_obj) & all(unlist(lapply(db_obj, function(x){inherits(x, "HibDatabase")}))) ), 
                          # Check model data
                          is.data.frame(mod_df),
                          # Check model
                          assertthat::assert_that(class(lme4_obj) == "lmerMod" || class(lme4_obj) == "glmerMod"), 
                          # Check coefficient df and metrics
                          is.data.frame(fixef_df),
                          is.character(ci_type) || is.na(ci_type),
                          is.numeric(metr_vec),
                          is.numeric(n_spp) || ( is.list(n_spp) & all(unlist(lapply(n_spp, function(x){is.numeric(x)}))) ),
                          is.numeric(n_taxa) || ( is.list(n_taxa) & all(unlist(lapply(n_taxa, function(x){is.numeric(x)}))) ),
                          is.numeric(n_sites) || ( is.list(n_sites) & all(unlist(lapply(n_sites, function(x){is.numeric(x)}))) )
                          )
                        
                        # ------ #
                        message(crayon::green("\nInitialising HibModelLME4 object ...") )
                        
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
                        # if list, get from each...
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
                        
                        
                        # Calc n sites/species used in  model - need to calc spp with div metrics!!
                        # N_sites <- dplyr::n_distinct(mod_df$SSBS)
                        
                      
                        # if ("SSBS" %in% colnames(mod_df)){
                        #   self$processing_flags$N_sites <- dplyr::n_distinct(mod_df$SSBS)
                        # } else if ("SSBS1" %in% colnames(mod_df)){
                        #   # print("SSBS1 found")
                        #   self$processing_flags$N_sites <- length(unique(c(mod_df$SSBS1, mod_df$SSBS2)))
                        # }
                        
                        # self$processing_flags$N_spp <- n_spp
                        
                        if (inherits(db_obj, "HibDatabase")){
                          self$dataset_name <- db_obj$dataset_name  
                          loc_df <- db_obj$sites[, c("SSBS", "Longitude" ,"Latitude")]

                        }
                        if (is.list(db_obj)){
                          self$dataset_name <- lapply(db_obj, function(x){x$dataset_name})
                          loc_df <- lapply(db_obj, function(x){x$sites[, c("SSBS", "Longitude" ,"Latitude")]}) |>
                            dplyr::bind_rows()
                        }
                        

                        # store site locations
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
                        # N spp/staxa/ites...
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
                        # stats <- list(number_sites = self$processing_flags$N_sites,
                        #               # number_lu = length( unique(self$data$Predominant_land_use) ),
                        #               # number_lui = length(unique(self$data$Use_intensity)),
                        #               number_species = paste(self$processing_flags$N_spp, collapse = ", ") ) #,
                        #               # filters_applied = paste(self$processing_flags$filtering_steps, collapse = "; "),
                        #               # div_metrics = paste(sort(unique(self$diversity$Metric)), collapse = ", ") )  
    
                        
                        message(
                          paste0("Hibiiscus fitted model (lme4)\n",
                                 "Last updated: ", crayon::green( self$last_update ),
                                 "\n",
                                 # TODO: Dataset name/s ?
                                 "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                                 # "\t Land use (class x intensity):  ", crayon::green(stats$number_lu), " x ",
                                 # crayon::green(stats$number_lui), "\n",
                                 "\t Number of species:             ", crayon::green(stats$number_species), "\n",
                                 "\t Number of taxa (orig. names):  ", crayon::green(stats$number_taxa), "\n",
                                 # n_taxa
                                 "\n",
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
                      
                      # #' @description
                      # #' Function to save the object to the specified file-path 
                      # #' ('output_path', must be '.rds').
                      # #'
                      # #' @param output_path (`character()`)
                      # #' @return (`NULL`).
                      # save_hib = function(output_path = NULL){
                      #   assertthat::assert_that(
                      #     is.null(output_path) || is.character(output_path)
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
                        
                      #   message(crayon::green(paste0("Saving hib model to: '", output_path, "'\n")) )
                        
                      #   saveRDS(self, output_path)
                      #   assertthat::assert_that(file.exists(output_path),
                      #                           msg = "Writing a copy of the model failed...")
                      #   invisible()
                      # },
                      
                      #' @description
                      #' Function to generate figures summarising the model.
                      #' 
                      #' @param type (`character()`) Options are: "fixed-effects",
                      #' "random-effects", "map-sites".
                      #' @return (`ggplot`).
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
                          # TODO: convert chr col of Coefficient_name to factor, based on order of terms in df
                          fixef_df$Coefficient_name <- factor(fixef_df$Coefficient_name,
                                                              levels = fixef_df$Coefficient_name)
  
                          # fixef_df$Coefficient_name <- factor(fixef_df$Coefficient_name)
                          # if ("(Intercept)" %in% fixef_df$Coefficient_name){
                          #   fixef_df$Coefficient <- relevel(fixef_df$Coefficient_name,
                          #                                   ref = "(Intercept)")
                          # }
                          
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
                          
                          rand_df <- as.data.frame(lme4::ranef(self$model)) |>
                            dplyr::mutate(lwr = condval - 1.96*condsd, 
                                          upr = condval + 1.96*condsd) |>
                            # order from small to large?
                            dplyr::group_by(grpvar, term) |>
                            dplyr::arrange(condval, .by_group = T) |>
                            dplyr::mutate(ID_ord = 1:dplyr::n()) |>
                            dplyr::ungroup() |>
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
                                             size = 0.1,
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
                      #' Function to assess the model.
                      #' 
                      #' @return (`NULL`).
                      check_model = function(){
                        m <- self$model
                        performance::check_model(m)
                      }
                    ),
                    private = list()
)