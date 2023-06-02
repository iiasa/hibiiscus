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
                            inherits(mod_obj, "hib_lme4_model") || inherits(mod_obj, "hib_inla_model") || inherits(mod_obj, "hib_brms_model"),
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
                                                          "fixedef_coefficients",
                                                          "extra_args"))
                            self$model_type <- "hib_brms_model"
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
                          if (self$model_type == "hib_lme4_model"){
                            m_form <- format(stats::formula(self$original_model_obj$model))
                          }
                          if (self$model_type == "hib_brms_model"){
                            m_form <- format(self$original_model_obj$model$formula$formula)
                          }
                          
                          
                          message(
                            paste0("Cross-validated ", self$model_type, " \n",
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
                          insamp_df <- self$insample_results$insample_metrics %>%
                            dplyr::filter(Metric_shrt %in% c("rmse", "rsq"))
                          
                          metr_tmp <- c("rmse" = insamp_df$Metric_val[insamp_df$Metric_shrt == "rmse"],
                                        "rsq"   = insamp_df$Metric_val[insamp_df$Metric_shrt == "rsq"])
                          
                          metr_tmp <- pad_fn(metr_tmp, 3, F)
                          names(metr_tmp) <- c("rmse", "rsq")
                          
                          message(
                            paste0("--------------------\n",
                                   "Cross-validation type:    ", crayon::green("In-sample"), "\n", 
                                   "No. of folds:             ", crayon::green("1"), "\n\n", 
                                   "RMSE:                     ", crayon::green( paste0( metr_tmp["rmse"], 
                                                                                        " (", 
                                                                                        round(insamp_df$Metric_lo[insamp_df$Metric_shrt == "rmse"],
                                                                                              3), 
                                                                                        ", ", 
                                                                                        round(insamp_df$Metric_hi[insamp_df$Metric_shrt == "rmse"],
                                                                                              3),  
                                                                                        ")" ) ) , "\n",
                                   "Predictive R-squared:     ", crayon::green( metr_tmp["rsq"] ) , "\n")
                          )
                          
                          
                          # per cv type
                          for (x in self$cv_results){
                            fld_typ <- x$fold_type
                            if (fld_typ == "Custom"){
                              fld_typ <- paste0(fld_typ, " (", x$fold_col, ")")
                            }
                            n_fld <- x$n_folds
                            
                            cv_df <- x$cv_metrics %>%
                              dplyr::filter(Metric_shrt %in% c("rmse", "rsq"))
                            
                            metr_tmp <- c("rmse" = cv_df$Metric_val[cv_df$Metric_shrt == "rmse"],
                                          "rsq"   = cv_df$Metric_val[cv_df$Metric_shrt == "rsq"])
                            
                            
                            metr_tmp <- pad_fn(metr_tmp, 3, F)
                            names(metr_tmp) <- c("rmse", "rsq")
                            
                            
                            message(
                              paste0("--------------------\n\n",
                                     "Cross-validation type:    ", crayon::green(fld_typ), "\n", 
                                     "No. of folds:             ", crayon::green(n_fld), "\n\n", 
                                     "RMSE:                     ", crayon::green(paste0( metr_tmp["rmse"], 
                                                                                         " (", 
                                                                                         round(cv_df$Metric_lo[cv_df$Metric_shrt == "rmse"],
                                                                                               3), 
                                                                                         ", ", 
                                                                                         round(cv_df$Metric_hi[cv_df$Metric_shrt == "rmse"],
                                                                                               3),  
                                                                                         ")" )), "\n",
                                     "Predictive R-squared:     ", crayon::green(metr_tmp["rsq"]) , "\n")
                            )
                          }
                          
                          
                          # maybe some original model metrics..
                          
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
                        
                        # plot function - 
                        plot = function(type = "metrics",
                                        metrics = c("rmse", "rsq")){
                          
                          # stop("Plotting functionality still in development...")
                          # Types:
                          #   -coefs - insample v cv (all folds v summary)
                          #   -metrics - over all fold types ; 
                          #              if multiple metrics, facet_wrap
                          #   -obs v pred, facet by fold_type, colour by fold_id
                          
                          assertthat::assert_that(is.character(type))
                          # Plot options
                          type <- match.arg(type,
                                            c("coefficients",
                                              "coefficients-summary",
                                              "obs-v-pred",
                                              "metrics"
                                            ))
                          
                          if (type == "metrics"){
                            metrics <- match.arg(metrics,
                                                 c("cor", "corsq",
                                                   "me", "mae", "rae", "mape", 
                                                   "smape", 
                                                   "mse", "rmse", "nrmse", 
                                                   "rmsle", "rse", "rrmse",
                                                   "rsq"),
                                                 several.ok = T)
                          }
                          
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
                          
                          if (type == "coefficients"){
                            # Get original and cv coefs from
                            
                            orig_coef <- self$original_model_obj$fixedef_coefficients
                            
                            cv_coef <- purrr::map_dfr(.x = self$cv_results,
                                                      .f = function(x){
                                                        tmp_df <- x$cv_coef_df
                                                        
                                                        cv_type <- x$fold_type
                                                        if (cv_type == "Custom"){
                                                          cv_type <- paste0(cv_type, " (", x$fold_col, ")")
                                                        }
                                                        
                                                        tmp_df$cv_type <- cv_type
                                                        
                                                        tmp_df$fold_id <- as.numeric(as.factor(tmp_df$fold_id))
                                                        
                                                        tmp_df
                                                      })
                            
                            
                            orig_coef <- orig_coef %>%
                              dplyr::mutate(coef_nm_num = as.numeric(as.factor(Coefficient_name))) %>%
                              dplyr::arrange(coef_nm_num)
                            
                            cv_coef <- cv_coef %>%
                              dplyr::mutate(coef_nm_num = as.numeric(as.factor(Coefficient_name)))
                            
                            # print(head(orig_coef))
                            
                            plt <- ggplot2::ggplot(cv_coef) +
                              ggplot2::geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                              
                              # %%%%
                              ggplot2::geom_rect(ggplot2::aes(xmin = coef_nm_num - 0.35,
                                                              xmax = coef_nm_num + 0.35,
                                                              ymin = loCI,
                                                              ymax = hiCI),
                                                 alpha = 0.3,
                                                 data = orig_coef) +
                              
                              ggplot2::geom_segment(ggplot2::aes(x = coef_nm_num - 0.35,
                                                                 xend = coef_nm_num + 0.35,
                                                                 y = Coefficient_value,
                                                                 yend = Coefficient_value),
                                                    alpha = 0.7,
                                                    data = orig_coef) +
                              # %%%%
                              
                              ggplot2::geom_errorbar(ggplot2::aes(x = coef_nm_num,
                                                                  ymin = loCI,
                                                                  ymax = hiCI,
                                                                  colour = as.character(fold_id)),
                                                     position = ggplot2::position_dodge(0.7),
                                                     width = 0,
                                                     alpha = 0.8) +
                              ggplot2::geom_point(ggplot2::aes(x = coef_nm_num,
                                                               y = Coefficient_value,
                                                               colour = as.character(fold_id)),
                                                  position = ggplot2::position_dodge(0.7),
                                                  alpha = 0.8) +
                              
                              # ggplot2::geom_errorbar(ggplot2::aes(x = Coefficient_name,
                              #                                     ymin = loCI,
                              #                                     ymax = hiCI),
                              #                        width = 0,
                              #                        alpha = 0.7,
                              #                        data = orig_coef) +
                              # ggplot2::geom_point(ggplot2::aes(x = Coefficient_name,
                              #                                  y = Coefficient_value),
                              #                     alpha = 0.7,
                              #                     data = orig_coef) +
                            
                            
                            ggplot2::labs(x = "Coefficient name",
                                          y = "Coefficient value") +
                              ggplot2::scale_x_continuous(breaks = orig_coef$coef_nm_num, 
                                                          labels = orig_coef$Coefficient_name) +
                              ggplot2::coord_flip() +
                              ggplot2::theme_bw() +
                              basic_thm +
                              ggplot2::facet_wrap(~cv_type) +
                              ggplot2::theme(legend.position = "none")
                            
                          }
                          
                          if (type == "coefficients-summary" ){
                            
                            orig_coef <- self$original_model_obj$fixedef_coefficients %>%
                              dplyr::mutate(cv_type = "In-sample",
                                            fold_id = NA)
                            
                            # summarise cv coefs to show range (95% quants?) 
                            cv_coef <- purrr::map_dfr(.x = self$cv_results,
                                                      .f = function(x){
                                                        tmp_df <- x$cv_coef_df
                                                        
                                                        cv_type <- x$fold_type
                                                        if (cv_type == "Custom"){
                                                          cv_type <- paste0(cv_type, " (", x$fold_col, ")")
                                                        }
                                                        
                                                        tmp_df$cv_type <- cv_type
                                                        
                                                        tmp_df$fold_id <- as.numeric(as.factor(tmp_df$fold_id))
                                                        
                                                        
                                                        tmp_df <- tmp_df %>%
                                                          dplyr::group_by(Coefficient_name, cv_type) %>%
                                                          dplyr::summarise(Coefficient_median = median(Coefficient_value),
                                                                           Coefficient_min = min(Coefficient_value),
                                                                           Coefficient_max = max(Coefficient_value),
                                                                           
                                                                           loCI_min = min(loCI),
                                                                           hiCI_max = max(hiCI)
                                                          )
                                                        
                                                      })
                            
                            orig_coef <- orig_coef %>%
                              dplyr::mutate(coef_nm_num = as.numeric(as.factor(Coefficient_name))) %>%
                              dplyr::arrange(coef_nm_num)
                            
                            cv_coef <- cv_coef %>%
                              dplyr::ungroup() %>%
                              dplyr::mutate(coef_nm_num = as.numeric(as.factor(Coefficient_name)))
                            
                            # print(head(orig_coef))
                            # print(head(cv_coef))
                            
                            plt <- ggplot2::ggplot(cv_coef) +
                              ggplot2::geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                              
                              # %%%%
                              ggplot2::geom_rect(ggplot2::aes(xmin = coef_nm_num - 0.35,
                                                              xmax = coef_nm_num + 0.35,
                                                              ymin = loCI,
                                                              ymax = hiCI),
                                                 alpha = 0.3,
                                                 data = orig_coef) +
                              
                              ggplot2::geom_segment(ggplot2::aes(x = coef_nm_num - 0.35,
                                                                 xend = coef_nm_num + 0.35,
                                                                 y = Coefficient_value,
                                                                 yend = Coefficient_value),
                                                    alpha = 0.7,
                                                    data = orig_coef) +
                              # %%%%
                              
                              ggplot2::geom_errorbar(ggplot2::aes(x = coef_nm_num,
                                                                  ymin = loCI_min,
                                                                  ymax = hiCI_max,
                                                                  colour = cv_type),
                                                     position = ggplot2::position_dodge(0.7),
                                                     width = 0,
                                                     size = 1,
                                                     alpha = 0.7) +
                              
                              ggplot2::geom_errorbar(ggplot2::aes(x = coef_nm_num,
                                                                  ymin = Coefficient_min,
                                                                  ymax = Coefficient_max,
                                                                  colour = cv_type),
                                                     position = ggplot2::position_dodge(0.7),
                                                     width = 0,
                                                     size = 1.75,
                                                     alpha = 0.7) +
                              
                              ggplot2::geom_point(ggplot2::aes(x = coef_nm_num,
                                                               y = Coefficient_median,
                                                               group = cv_type),
                                                  colour = "white",
                                                  position = ggplot2::position_dodge(0.7)) +
                              
                              # ggplot2::geom_errorbar(ggplot2::aes(x = Coefficient_name,
                              #                                     ymin = loCI,
                              #                                     ymax = hiCI),
                              #                        # position = ggplot2::position_nudge(x = -0.325),
                              #                        width = 0,
                              #                        alpha = 0.7,
                              #                        data = orig_coef) +
                              # ggplot2::geom_point(ggplot2::aes(x = Coefficient_name,
                              #                                  y = Coefficient_value),
                              #                     # position = ggplot2::position_nudge(x = -0.325),
                            #                     alpha = 0.7,
                            #                     data = orig_coef) +
                            
                            ggplot2::scale_x_continuous(breaks = orig_coef$coef_nm_num, 
                                                        labels = orig_coef$Coefficient_name) +
                              ggplot2::labs(x = "Coefficient name",
                                            y = "Coefficient value",
                                            colour = "Cross-\nvalidation\ntype") +
                              ggplot2::coord_flip() +
                              ggplot2::theme_bw() +
                              basic_thm 
                          }
                          if (type == "obs-v-pred"){
                            # Get obs and pred - across cv + insample
                            insamp_df <- self$insample_results$insample_pred_df %>%
                              dplyr::mutate(cv_type = "In-sample",
                                            fold_id = NA)
                            
                            cv_df <- purrr::map_dfr(.x = self$cv_results,
                                                    .f = function(x){
                                                      tmp_df <- x$cv_pred_df
                                                      
                                                      cv_type <- x$fold_type
                                                      if (cv_type == "Custom"){
                                                        cv_type <- paste0(cv_type, " (", x$fold_col, ")")
                                                      }
                                                      
                                                      tmp_df$cv_type <- cv_type
                                                      
                                                      tmp_df$fold_id <- as.numeric(as.factor(tmp_df$fold_id))
                                                      
                                                      tmp_df
                                                    })
                            
                            
                            df <- dplyr::bind_rows(insamp_df,
                                                   cv_df) %>%
                              dplyr::mutate(cv_type = factor(cv_type),
                                            cv_type = relevel(cv_type, ref = "In-sample"))
                            
                            rng <- range(c(df$y, df$prediction), na.rm = T)
                            
                            
                            plt <- ggplot2::ggplot(df) +
                              ggplot2::geom_abline(intercept = 0, 
                                                   slope = 1,
                                                   lty = 2, colour = "darkgrey") +
                              ggplot2::geom_point(ggplot2::aes(x = y,
                                                               y = prediction,
                                                               colour = as.character(fold_id)),
                                                  alpha = 0.25) +
                              ggplot2::labs(x = "Observed",
                                            y = "Predicted") +
                              ggplot2::coord_equal(xlim = rng,
                                                   ylim = rng) +
                              ggplot2::theme_bw() +
                              basic_thm +
                              ggplot2::facet_wrap(~cv_type) +
                              ggplot2::theme(legend.position = "none")
                          }
                          
                          if (type == "metrics"){
                            # Match metrics vec to names in metrics data
                            
                            cv_df <- purrr::map_dfr(.x = self$cv_results,
                                                    .f = function(x){
                                                      tmp_df <- x$cv_metrics
                                                      
                                                      cv_type <- x$fold_type
                                                      if (cv_type == "Custom"){
                                                        cv_type <- paste0(cv_type, " (", x$fold_col, ")")
                                                      }
                                                      tmp_df$cv_type <- cv_type
                                                      
                                                      tmp_df
                                                    }) %>%
                              dplyr::filter(Metric_shrt %in% metrics) %>%
                              dplyr::mutate(cv_type_num = as.numeric(as.factor(cv_type)))
                            
                            insamp_df <- self$insample_results$insample_metrics %>%
                              dplyr::filter(Metric_shrt %in% metrics)
                            
                            # plot metrics by cv type, facet per metric, colour per cv type
                            plt <- ggplot2::ggplot(cv_df) +
                              
                              ggplot2::geom_segment(ggplot2::aes(x = min(cv_df$cv_type_num),
                                                                 xend = max(cv_df$cv_type_num),
                                                                 y = Metric_val,
                                                                 yend = Metric_val),
                                                    data = insamp_df) +
                              ggplot2::geom_rect(ggplot2::aes(xmin = min(cv_df$cv_type_num),
                                                              xmax = max(cv_df$cv_type_num),
                                                              ymin = Metric_lo,
                                                              ymax = Metric_hi),
                                                 alpha = 0.3,
                                                 data = insamp_df) +
                              
                              ggplot2::geom_errorbar(ggplot2::aes(x = cv_type_num,
                                                                  ymin = Metric_lo,
                                                                  ymax = Metric_hi,
                                                                  colour = cv_type),
                                                     width = 0) +
                              ggplot2::geom_point(ggplot2::aes(x = cv_type_num,
                                                               y = Metric_val,
                                                               colour = cv_type),
                                                  size = 2) +
                              ggplot2::labs(x = ggplot2::element_blank(),
                                            y = "Metric value",
                                            colour = "Cross-\nvalidation\ntype") +
                              ggplot2::theme_bw() +
                              basic_thm +
                              ggplot2::facet_wrap(~Metric_fct, scales = "free_y",
                                         labeller = ggplot2::label_parsed) +
                              ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                             axis.ticks.x = ggplot2::element_blank(),
                                             aspect.ratio = 1)
                            
                            # define metr_levs and metr_labs once...
                            metr_levs <- c("Correlation",
                                           "Correlation_squared",
                                           
                                           "Mean_error", 
                                           "Mean_absolute_error", 
                                           "Relative_absolute_error", 
                                           "Mean_absolute_percentage_error", 
                                           "Symmetric_mean_absolute_percentage_error", 
                                           
                                           "Mean_squared_error", 
                                           "Root_mean_squared_error", 
                                           "Normalised_root_mean_squared_error",
                                           "Relative_root_mean_squared_error",
                                           "Root_mean_squared_log_error",
                                           "Relative_squared_error", 
                                           
                                           "Predictive_R_squared")
                            
                            metr_labs <- c("Correlation",
                                           "Correlation^2",
                                           
                                           "ME~(Bias)", 
                                           "MAE", 
                                           "RAE", 
                                           "MAPE", 
                                           "sMAPE", 
                                           
                                           "MSE", 
                                           "RMSE", 
                                           "NRMSE",
                                           "RRMSE",
                                           "RMSLE",
                                           "RSE", 
                                           
                                           "Predictive~R^2")
                            
                            # if cor/cor2/rs2, add lines at 0, bad, 1, good
                            if (any(c("cor", "corsq", "rsq") %in% metrics)){
                              
                              cor_df <- data.frame(bad = 0,
                                                   good = 1,
                                                   Metric = c("Correlation",
                                                              "Correlation_squared",
                                                              "Predictive_R_squared"))
                              
                              cor_df <- cor_df %>%
                                dplyr::filter(Metric %in% cv_df$Metric) %>%
                                dplyr::mutate(Metric_fct = factor(Metric,
                                                                  levels = metr_levs,
                                                                  labels = metr_labs ))
                              
                              plt <- plt +
                                ggplot2::geom_hline(ggplot2::aes(yintercept = bad), 
                                                    lty = 2, 
                                                    colour = "darkgrey",
                                                    data = cor_df) +
                                ggplot2::geom_hline(ggplot2::aes(yintercept = good), 
                                                    lty = 1, 
                                                    colour = "darkgrey",
                                                    data = cor_df) 
                            }
                            # for all others, add line at 0, good
                            if (any(c("me", "mae", "rae", "mape", 
                                      "smape", 
                                      "mse", "rmse", "nrmse", 
                                      "rrmse", "rmsle", "rse") %in% metrics)){
                              err_df <- data.frame(good = 0,
                                                   Metric = c("Mean_error", 
                                                              "Mean_absolute_error", 
                                                              "Relative_absolute_error", 
                                                              "Mean_absolute_percentage_error", 
                                                              "Symmetric_mean_absolute_percentage_error", 
                                                              
                                                              "Mean_squared_error", 
                                                              "Root_mean_squared_error", 
                                                              "Normalised_root_mean_squared_error",
                                                              "Relative_root_mean_squared_error",
                                                              "Root_mean_squared_log_error",
                                                              "Relative_squared_error"))
                              
                              err_df <- err_df %>%
                                dplyr::filter(Metric %in% cv_df$Metric) %>%
                                dplyr::mutate(Metric_fct = factor(Metric,
                                                                  levels = metr_levs,
                                                                  labels = metr_labs ))
                              
                              plt <- plt +
                                ggplot2::geom_hline(ggplot2::aes(yintercept = good), 
                                                    lty = 1, 
                                                    colour = "darkgrey",
                                                    data = err_df) 
                            }
                            
                            if ("rse" %in% metrics){
                              rse_df <- data.frame(bad = 1,
                                                   Metric = c("Relative_squared_error"))
                              
                              rse_df <- rse_df %>%
                                dplyr::mutate(Metric_fct = factor(Metric,
                                                                  levels = metr_levs,
                                                                  labels = metr_labs ))
                              
                              plt <- plt +
                                ggplot2::geom_hline(ggplot2::aes(yintercept = bad), 
                                                    lty = 2, 
                                                    colour = "darkgrey",
                                                    data = rse_df) 
                            }
                            
                          }
                          
                          (plt)
                        }
                      ),
                      private = list()
)
