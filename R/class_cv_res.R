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
                                          "model", "processing_flags",
                                          "fixedef_coefficients")
                          if (inherits(mod_obj, "hib_lme4_model")){ 
                            core_items <- c(core_items, c("metrics", 
                                                          "ci_type"
                                                          # "fixedef_coefficients"
                                                          ))
                            self$model_type <- "hib_lme4_model"
                          }
                          if (inherits(mod_obj, "hib_inla_model")){ 
                            core_items <- c(core_items, c("user_formula",
                                                          "inla_fitted_formula", 
                                                          "inla_incl_spatial_term",
                                                          # "inla_mesh",
                                                          "inla_mesh_settings"))
                            self$model_type <- "hib_inla_model"
                          }
                          if (inherits(mod_obj, "hib_brms_model")){ 
                            core_items <- c(core_items, c("metrics", 
                                                          # "fixedef_coefficients",
                                                          "extra_args"))
                            self$model_type <- "hib_brms_model"
                          }
                          
                          
                          self$original_model_obj <- purrr::map(.x = core_items, 
                                                                .f = function(x, 
                                                                              obj = mod_obj){
                                                                  obj[[x]]
                                                                })
                          
                          names(self$original_model_obj) <- core_items
                          
                          # If inla, add fixedef_coefficients
                          # TODO: add this to hib_inla_model class - and rm below
                          # if (inherits(mod_obj, "hib_inla_model")){
                          #   fixef_df <- mod_obj$model$summary.fixed %>%
                          #     tibble::rownames_to_column("Coefficient_name") %>%
                          #     dplyr::rename(Coefficient_value = mean,
                          #                   loCI = `0.025quant`,
                          #                   median = `0.5quant`,
                          #                   hiCI = `0.975quant`)
                          #   self$original_model_obj$fixedef_coefficients <- fixef_df
                          # }
                          
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
                          if (self$model_type == "hib_inla_model"){
                            # m_form <- format(self$original_model_obj$user_formula$full)
                            # m_form <- format(self$original_model_obj$model$.args$formula)
                            m_form <- as.character(self$original_model_obj$user_formula$full)
                            m_form <- paste(form_chr[2], form_chr[1], form_chr[3])
                            if (self$original_model_obj$inla_incl_spatial_term){
                              m_form <- paste(m_form,
                                              "+ f(s, model = site.spde)")
                            }
                          }
                          
                          if (self$model_type == "hib_brms_model"){
                            m_form <- format(self$original_model_obj$model$formula$formula)
                          }
                          
                          # TODO: additional metrics: cor2, me (bias)
                          df1 <- self$insample_results$insample_metrics %>%
                            dplyr::filter(Metric_shrt %in% c("rmse", "rsq")) 
                          
                          insamp_df <- data.frame(CV_type = c("In-sample"),
                                                  N_folds = c(1),
                                                  RMSE    = c(df1$Metric_val[df1$Metric_shrt == "rmse"]),
                                                  RMSE_lo = c(df1$Metric_lo[df1$Metric_shrt == "rmse"]),
                                                  RMSE_hi = c(df1$Metric_hi[df1$Metric_shrt == "rmse"]),
                                                  Predictive_R2 = c(df1$Metric_val[df1$Metric_shrt == "rsq"]))
                          
                          
                          cv_df <- purrr::map_dfr(self$cv_results,
                                                  function(x){
                                                    fld_typ <- x$fold_type
                                                    if (fld_typ == "Custom"){
                                                      fld_typ <- paste0(fld_typ, " (", x$fold_col, ")")
                                                    }
                                                    n_fld <- x$n_folds
                                                    
                                                    cv_df <- x$cv_metrics %>%
                                                      dplyr::filter(Metric_shrt %in% c("rmse", "rsq"))
                                                    
                                                    out_df <- data.frame(CV_type = c(fld_typ),
                                                                         N_folds = c(n_fld),
                                                                         RMSE    = c(cv_df$Metric_val[cv_df$Metric_shrt == "rmse"]),
                                                                         RMSE_lo = c(cv_df$Metric_lo[cv_df$Metric_shrt == "rmse"]),
                                                                         RMSE_hi = c(cv_df$Metric_hi[cv_df$Metric_shrt == "rmse"]),
                                                                         Predictive_R2 = c(cv_df$Metric_val[cv_df$Metric_shrt == "rsq"]))
                                                    
                                                  })
                          
                          
                          df <- dplyr::bind_rows(insamp_df,
                                                 cv_df)
                          # drop n from cv_type for random
                          df$CV_type[grepl("Random", df$CV_type)] <- "Random"
                          df$RMSE <- pad_fn(df$RMSE, 3, T)
                          df$RMSE_lo <- pad_fn(df$RMSE_lo, 3, T)
                          df$RMSE_hi <- pad_fn(df$RMSE_hi, 3, T)
                          df$Predictive_R2 <- pad_fn(df$Predictive_R2, 3, T)
                          
                          
                          df$RMSE <- paste0(df$RMSE, " (", df$RMSE_lo, ", ", df$RMSE_hi, ")")
                          
                          
                          message(
                            paste0("Cross-validated ", self$model_type, " \n\n",
                                   "(Last updated: ", crayon::green( self$last_update ), ")",
                                   # "\n",
                                   # "\t Cross-validation types:        ", crayon::green(stats$number_sites), "\n"
                                   # "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                                   # # "\t Land use (class x intensity):  ", crayon::green(stats$number_lu), " x ",
                                   # # crayon::green(stats$number_lui), "\n",
                                   # "\t Number of species:             ", crayon::green(stats$number_species), "\n",
                                   "\n\n",
                                   "--------------------\n\n",
                                   "Summary of cross-validated model\n\n",
                                   # "Formula:                  ", crayon::green( m_form ) , "\n"
                                   "Formula:\n", 
                                   crayon::green( m_form ) , "\n\n",
                                   "Performance metrics:\n"
                            )
                          )
                          
                          
                          df %>%
                            dplyr::select(CV_type, N_folds, RMSE, Predictive_R2) %>%
                            dplyr::rename(`CV type` = CV_type,
                                          `N folds` = N_folds,
                                          `Predictive R2` = Predictive_R2) %>%
                            insight::format_table(use_symbols = T, digits = 0) %>%
                            insight::export_table() %>%
                            message()
                          
                          
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
                                              "obs-v-error", 
                                              "pred-v-error",
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
                          
                          # TODO
                          # specify order of cv_type levels 
                          base_cv_levs <- c("In-sample", "Random-2", "Random-3", 
                                            "Random-5", "Random-10",
                                            "Source_ID", "SS", #"SSBS",
                                            "Country", "UN_subregion", "UN_region",
                                            "Biome", "Realm")
                          # plus add unique custom cv types...
                          
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
                            
                            # Add cv_type that isn't in base_cv_levs to cv_levs
                            add_cv_levs <- unique(cv_coef$cv_type)
                            add_cv_levs <- add_cv_levs[!(add_cv_levs %in% base_cv_levs)]
                            cv_levs <- c(base_cv_levs, add_cv_levs)
                            # set levs for cv_type
                            
                            orig_coef <- orig_coef %>%
                              dplyr::mutate(coef_nm_num = as.numeric(as.factor(Coefficient_name))) %>%
                              dplyr::arrange(coef_nm_num)
                            
                            cv_coef <- cv_coef %>%
                              dplyr::mutate(coef_nm_num = as.numeric(as.factor(Coefficient_name)),
                                            cv_type = factor(cv_type, 
                                                             levels = cv_levs))
                            
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
                                                     alpha = 0.7) +
                              ggplot2::geom_point(ggplot2::aes(x = coef_nm_num,
                                                               y = Coefficient_value,
                                                               colour = as.character(fold_id)),
                                                  position = ggplot2::position_dodge(0.7),
                                                  alpha = 0.7) +
                              
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
                            
                            # Add cv_type that isn't in base_cv_levs to cv_levs
                            add_cv_levs <- unique(cv_coef$cv_type)
                            add_cv_levs <- add_cv_levs[!(add_cv_levs %in% base_cv_levs)]
                            cv_levs <- c(base_cv_levs, add_cv_levs)
                            
                            
                            orig_coef <- orig_coef %>%
                              dplyr::mutate(coef_nm_num = as.numeric(as.factor(Coefficient_name))) %>%
                              dplyr::arrange(coef_nm_num)
                            
                            cv_coef <- cv_coef %>%
                              dplyr::ungroup() %>%
                              dplyr::mutate(coef_nm_num = as.numeric(as.factor(Coefficient_name)),
                                            cv_type = factor(cv_type, 
                                                             levels = cv_levs))
                            
                            # print(head(orig_coef))
                            # print(head(cv_coef))
                            
                            plt <- ggplot2::ggplot(cv_coef) +
                              ggplot2::geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                              
                              # %%%%
                              ggplot2::geom_rect(ggplot2::aes(xmin = coef_nm_num - 0.35,
                                                              xmax = coef_nm_num + 0.35,
                                                              ymin = loCI,
                                                              ymax = hiCI),
                                                 alpha = 0.35,
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
                                                     size = 0.9,
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
                          
                          # TODO: 
                          if (type %in% c("obs-v-pred", "obs-v-error", "pred-v-error") ){
                            
                            # Get obs and pred and error - across cv + insample
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
                            
                            
                            # Add cv_type that isn't in base_cv_levs to cv_levs
                            add_cv_levs <- unique(cv_df$cv_type)
                            add_cv_levs <- add_cv_levs[!(add_cv_levs %in% base_cv_levs)]
                            cv_levs <- c(base_cv_levs, add_cv_levs)
                            
                            df <- dplyr::bind_rows(insamp_df,
                                                   cv_df) %>%
                              dplyr::mutate(cv_type = factor(cv_type,
                                                             levels = cv_levs))
                            
                            
                            if (type == "obs-v-pred"){
                              rng <- range(c(df$y, df$prediction), na.rm = T)  
                            } else {
                              y_abs_max <- max(abs(range(df$error, na.rm = T)))
                            }
                            
                            
                            plt <- ggplot2::ggplot(df)
                            
                            if (type == "obs-v-pred"){
                              plt <- plt +
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
                                                     ylim = rng) 
                            } else {
                              plt <- plt +
                                ggplot2::geom_hline(yintercept = 0, 
                                                     lty = 2, colour = "darkgrey") +
                                ggplot2::scale_y_continuous(limits = c(-y_abs_max, y_abs_max))
                            }
                            
                            if (type == "obs-v-error"){
                              plt <- plt +
                                ggplot2::geom_point(ggplot2::aes(x = y,
                                                                 y = error,
                                                                 colour = as.character(fold_id)),
                                                    alpha = 0.25) +
                                ggplot2::labs(x = "Observed",
                                              y = "Error")
                            }
                            if (type == "pred-v-error"){
                              plt <- plt +
                                ggplot2::geom_point(ggplot2::aes(x = prediction,
                                                                 y = error,
                                                                 colour = as.character(fold_id)),
                                                    alpha = 0.25) +
                                ggplot2::labs(x = "Predicted",
                                              y = "Error")
                            }
                            plt <- plt +
                              ggplot2::theme_bw() +
                              basic_thm +
                              ggplot2::facet_wrap(~cv_type) +
                              ggplot2::theme(legend.position = "none")
                            
                          }
                          
                          if (type == "metrics"){
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
                              dplyr::filter(Metric_shrt %in% metrics)
                            
                            
                            # %%%%
                            cv_fold_df <- purrr::map_dfr(.x = self$cv_results,
                                                    .f = function(x){
                                                      tmp_df <- x$cv_fold_metrics
                                                      
                                                      cv_type <- x$fold_type
                                                      if (cv_type == "Custom"){
                                                        cv_type <- paste0(cv_type, " (", x$fold_col, ")")
                                                      }
                                                      tmp_df$cv_type <- cv_type
                                                      
                                                      tmp_df$fold_id <- as.numeric(tmp_df$fold_id)
                                                      
                                                      tmp_df
                                                    }) %>%
                              dplyr::filter(Metric_shrt %in% metrics)
                            # %%%%
                            
                            insamp_df <- self$insample_results$insample_metrics %>%
                              dplyr::filter(Metric_shrt %in% metrics)
                            
                            # print(insamp_df)
                            # Add cv_type that isn't in base_cv_levs to cv_levs
                            add_cv_levs <- unique(cv_df$cv_type)
                            add_cv_levs <- add_cv_levs[!(add_cv_levs %in% base_cv_levs)]
                            
                            cv_levs <- base_cv_levs[base_cv_levs %in% unique(cv_df$cv_type)]
                            cv_levs <- c(cv_levs, add_cv_levs)
                            
                            cv_df <- cv_df %>%
                              dplyr::mutate(cv_type = factor(cv_type,
                                                             levels = cv_levs),
                                            cv_type_num = as.numeric(cv_type))
                            
                            cv_fold_df <- cv_fold_df %>%
                              dplyr::mutate(cv_type = factor(cv_type,
                                                             levels = cv_levs),
                                            cv_type_num = as.numeric(cv_type))
                            # print(cv_df)
                            # plot metrics by cv type, facet per metric, colour per cv type
                            plt <- ggplot2::ggplot(cv_df)
                            
                            
                            # %%%%
                            # Add some baselines
                            # %%%%
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
                            
                            plt <- plt +  
                              ggplot2::geom_segment(ggplot2::aes(x = min(cv_df$cv_type_num)-0.5,
                                                                 xend = max(cv_df$cv_type_num)+0.5,
                                                                 y = Metric_val,
                                                                 yend = Metric_val),
                                                    data = insamp_df) +
                              ggplot2::geom_rect(ggplot2::aes(xmin = min(cv_df$cv_type_num)-0.5,
                                                              xmax = max(cv_df$cv_type_num)+0.5,
                                                              ymin = Metric_lo,
                                                              ymax = Metric_hi),
                                                 alpha = 0.3,
                                                 data = insamp_df) +
                              
                              ggplot2::geom_point(ggplot2::aes(x = cv_type_num,
                                                               y = Metric_val),
                                                  # , colour = cv_type
                                                  alpha = 0.5,
                                                  size = 0.75,
                                                  position = ggplot2::position_jitter(width = 0.2),
                                                  data = cv_fold_df) +
                              
                              ggplot2::geom_errorbar(ggplot2::aes(x = cv_type_num,
                                                                  ymin = Metric_lo,
                                                                  ymax = Metric_hi,
                                                                  colour = cv_type),
                                                     width = 0) +
                              ggplot2::geom_point(ggplot2::aes(x = cv_type_num,
                                                               y = Metric_val,
                                                               colour = cv_type),
                                                  size = 2.5,
                                                  alpha = 0.75) +
                              
                              
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
                            
                            
                            
                          }
                          
                          (plt)
                        }
                      ),
                      private = list()
)
