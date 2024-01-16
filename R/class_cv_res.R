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
HibCV <- R6Class("HibCV",
                      # "HibCV"
                      # TODO: @export

                      # Public data objects and functions
                      public = list(

                        #' @field last_update (`character`)\cr
                        #' Identifies the most recent update to the object (YYYY-MM-DD).
                        last_update = NULL,
                        
                        #' @field model_type (`character`)\cr
                        #' Identifies the type of model; one of "HibModelLME4", 
                        #' "HibModelINLA", "HibModelBRMS".
                        model_type = NULL,
                        
                        #' @field model_type (`list`)\cr
                        #' Contains a list of core model information
                        original_model_obj = NULL, # All model info, except data, and functions
                        
                        #' @field cv_results (`list`)\cr
                        #' Contains the list of results from cross-validation
                        cv_results = NULL,
                        
                        #' @field insample_results (`list`)\cr
                        #' Contains the list of results from in-sample validation
                        insample_results = NULL,
                        
                        #' @field re_form (`list`)\cr
                        #' A list of the re_form model structures used in evaluation.
                        re_form = NULL,
                        
                        #' @field use_original_mesh (`logical` or `NULL`)\cr
                        #' Indicates use of original mesh for spatial INLA models.
                        use_original_mesh = NULL,
                        
                        # Executed when created
                        #' @description
                        #' Creates a new instance of this [R6][R6::R6Class] class.
                        #' Note, this object is constructed via the [hib_cv_model]
                        #' function found within this package.
                        #' @param mod_obj (`HibModelXXXX`)\cr
                        #'   The model object used in cross-validation.
                        #' @param cv_res (`list`)\cr
                        #'   The results of cross-validation.
                        #' @param insamp_res (`list`)\cr
                        #'   The results of ins-sample validation.
                        #' @param re_form (`list`)\cr
                        #'   The list of random effects structures used in model evaluation.
                        #' @param use_orig_mesh (`logical` or `NULL`)\cr
                        #'   Indicates use of original mesh for spatial INLA models
                        initialize = function(mod_obj,
                                              cv_res,
                                              insamp_res,
                                              re_form,
                                              use_orig_mesh){
                          # - some assertion checks
                          assertthat::assert_that(
                            # Check model
                            inherits(mod_obj, "HibModelLME4") || 
                              inherits(mod_obj, "HibModelINLA") || 
                              inherits(mod_obj, "HibModelBRMS"),
                            # Check cv outputs
                            is.list(cv_res),
                            is.list(insamp_res),
                            is.list(re_form),
                            is.logical(use_orig_mesh) || is.null(use_orig_mesh)
                          )
                          
                          # ------ #
                          message(crayon::green("\nInitialising HibCV object ...") )
                          
                          # %%%%%%%%%%%%%%%
                          
                          # Save inputs in self
                          # dataset_name
                          core_items <- c("last_update", "dataset_name",
                                          "model", "processing_flags",
                                          "fixedef_coefficients",
                                          "site_locations")
                          # 
                          
                          if (inherits(mod_obj, "HibModelLME4")){ 
                            core_items <- c(core_items, c("metrics", 
                                                          "ci_type"
                                                          # "fixedef_coefficients"
                                                          ))
                            self$model_type <- "HibModelLME4"
                          }
                          if (inherits(mod_obj, "HibModelINLA")){ 
                            core_items <- c(core_items, c("user_formula",
                                                          "inla_fitted_formula", 
                                                          "inla_incl_spatial_term",
                                                          # "inla_mesh",
                                                          "inla_mesh_settings",
                                                          "inla_mesh_source"))
                            self$model_type <- "HibModelINLA"
                          }
                          if (inherits(mod_obj, "HibModelBRMS")){ 
                            core_items <- c(core_items, c("metrics", 
                                                          # "fixedef_coefficients",
                                                          "extra_args"))
                            self$model_type <- "HibModelBRMS"
                          }
                          
                          
                          self$original_model_obj <- purrr::map(.x = core_items, 
                                                                .f = function(x, 
                                                                              obj = mod_obj){
                                                                  obj[[x]]
                                                                })
                          
                          names(self$original_model_obj) <- core_items
                          
                          
                          self$re_form <- re_form

                          if (is.null(use_orig_mesh)){
                            self$use_original_mesh <- NA
                          } else {
                            self$use_original_mesh <- use_orig_mesh
                          }
                          
                          # If inla, add fixedef_coefficients
                          # TODO: add this to hib_inla_model class - and rm below
                          # if (inherits(mod_obj, "hib_inla_model")){
                          #   fixef_df <- mod_obj$model$summary.fixed |>
                          #     tibble::rownames_to_column("Coefficient_name") |>
                          #     dplyr::rename(Coefficient_value = mean,
                          #                   loCI = `0.025quant`,
                          #                   median = `0.5quant`,
                          #                   hiCI = `0.975quant`)
                          #   self$original_model_obj$fixedef_coefficients <- fixef_df
                          # }
                          
                          self$cv_results <- cv_res
                          self$insample_results <- insamp_res
                          
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
                          
                          # stats <- list(number_sites = self$original_model_obj$processing_flags$N_sites,
                          #               number_species = paste(self$original_model_obj$processing_flags$N_spp, collapse = ", ") ) #,
                          
                          
                          
                          # Model form
                          if (self$model_type == "HibModelLME4"){
                            m_form <- format(stats::formula(self$original_model_obj$model))
                          }
                          if (self$model_type == "HibModelINLA"){
                            # m_form <- format(self$original_model_obj$user_formula$full)
                            # m_form <- format(self$original_model_obj$model$.args$formula)
                            m_form <- as.character(self$original_model_obj$user_formula$full)
                            # m_form <- paste(form_chr[2], form_chr[1], form_chr[3])
                            m_form <- paste(m_form[2], m_form[1], m_form[3])
                            if (self$original_model_obj$inla_incl_spatial_term){
                              m_form <- paste(m_form,
                                              "+ f(s, model = site.spde)")
                            }
                          }
                          
                          if (self$model_type == "HibModelBRMS"){
                            m_form <- format(self$original_model_obj$model$formula$formula)
                          }
                          
                          # TODO: additional metrics: cor2, me (bias)
                          # df1 <- self$insample_results$insample_metrics |>
                          #   dplyr::filter(Metric_shrt %in% c("rmse", "rsq")) 
                          
                          # TODO: check nrow of cv_results[[i]]$cv_metrics
                          # insample_results$insamp_metrics/


                          # TODO: Incl Random_effect_form column - better approach
                          insamp_df <- self$insample_results$insample_metrics |>
                            dplyr::filter(Metric_shrt %in% c("rmse")) |>
                            dplyr::rename(RMSE = Metric_val,
                                          RMSE_lo = Metric_lo,
                                          RMSE_hi = Metric_hi) |>
                            dplyr::select(RMSE, RMSE_lo, RMSE_hi, Random_effect_form, N) |>
                            dplyr::left_join(self$insample_results$insample_metrics |>
                                               dplyr::filter(Metric_shrt %in% c("rsq")) |>
                                               dplyr::rename(Predictive_R2 = Metric_val) |>
                                               dplyr::select(Predictive_R2, Random_effect_form),
                                             by = "Random_effect_form"
                            ) |>
                            dplyr::mutate(CV_type = "In-sample",
                                          N_folds = 1)
                          
                          
                          # cv_specs - expand.grid(re_from, fld_type)
                          # Any not in cv_df - add with NAs for metrics and N = 0
                          cv_specs <- purrr::map_dfr(self$cv_results,
                                                     function(x, re_f = as.character(self$re_form)){
                                                       fld_typ <- x$fold_type
                                                       if (fld_typ == "Custom"){
                                                         fld_typ <- paste0(fld_typ, " (", x$fold_col, ")")
                                                       }
                                                       n_fld <- x$n_folds
                                                       
                                                       return(data.frame(CV_type = fld_typ,
                                                                         N_folds = n_fld,
                                                                         Random_effect_form = re_f))
                                                     })
                          cv_specs$CV_type[grepl("Random", cv_specs$CV_type)] <- "Random"

                          cv_df <- purrr::map_dfr(self$cv_results,
                                                  function(x){
                                                    fld_typ <- x$fold_type
                                                    if (fld_typ == "Custom"){
                                                      fld_typ <- paste0(fld_typ, " (", x$fold_col, ")")
                                                    }
                                                    n_fld <- x$n_folds
                                                    
                                                    if (is.null(x$cv_metrics)){
                                                      return(NULL)
                                                      # TODO: return and df with unavailable fold_type and re_form
                                                    }
                                                    if (n_fld < 2){
                                                      return(NULL)
                                                    }
                                                    # cv_df <- x$cv_metrics |>
                                                    #   dplyr::filter(Metric_shrt %in% c("rmse", "rsq"))
                                                    
                                                    # out_df <- data.frame(CV_type = c(fld_typ),
                                                    #                      N_folds = c(n_fld),
                                                    #                      Random_effect_form = c(cv_df$Random_effect_form[cv_df$Metric_shrt == "rmse"]),
                                                    #                      RMSE    = c(cv_df$Metric_val[cv_df$Metric_shrt == "rmse"]),
                                                    #                      RMSE_lo = c(cv_df$Metric_lo[cv_df$Metric_shrt == "rmse"]),
                                                    #                      RMSE_hi = c(cv_df$Metric_hi[cv_df$Metric_shrt == "rmse"]),
                                                    #                      Predictive_R2 = c(cv_df$Metric_val[cv_df$Metric_shrt == "rsq"]))
                                                    
                                                    out_df <- x$cv_metrics |>
                                                      dplyr::filter(Metric_shrt %in% c("rmse")) |>
                                                      dplyr::rename(RMSE = Metric_val,
                                                                    RMSE_lo = Metric_lo,
                                                                    RMSE_hi = Metric_hi) |>
                                                      dplyr::select(RMSE, RMSE_lo, RMSE_hi, Random_effect_form, N) |>
                                                      dplyr::left_join(x$cv_metrics |>
                                                                         dplyr::filter(Metric_shrt %in% c("rsq")) |>
                                                                         dplyr::rename(Predictive_R2 = Metric_val) |>
                                                                         dplyr::select(Predictive_R2, Random_effect_form),
                                                                       by = "Random_effect_form") |>
                                                      dplyr::mutate(CV_type = fld_typ,
                                                                    N_folds = n_fld)
                                                    # if (n_fld > 1){
                                                    #   
                                                    # } else {
                                                    #   out_df <- data.frame(CV_type = c(fld_typ),
                                                    #                        N_folds = c(n_fld),
                                                    #                        RMSE    = NA,
                                                    #                        RMSE_lo = NA,
                                                    #                        RMSE_hi = NA,
                                                    #                        Predictive_R2 = NA)
                                                    # }
                                                    
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
                          
                          missing_cv_specs <- dplyr::anti_join(cv_specs,
                                                               df, 
                                                               by = c("CV_type", "Random_effect_form"))
                          
                          if (nrow(missing_cv_specs) > 0){
                            missing_cv_specs <- missing_cv_specs |>
                              dplyr::mutate(RMSE = "N/A",
                                            Predictive_R2 = "N/A",
                                            N = 0)
                            df <- df |>
                              dplyr::bind_rows(missing_cv_specs)
                          }
                          
                          mtype_shrt <- NA 
                          if (self$model_type == "HibModelLME4"){
                            mtype_shrt <- "lme4"
                          }
                          if (self$model_type == "HibModelBRMS"){
                            mtype_shrt <- "brms"
                          }
                          if (self$model_type == "HibModelLME4"){
                            mtype_shrt <- "INLA"
                          }
                          message(
                            paste0("Hibiiscus cross-validated model (", mtype_shrt, ")\n",
                                   "Last updated: ", crayon::green( self$last_update ),
                                   "\n\n",
                                   "--------------------\n\n",
                                   "Summary of cross-validated model\n\n",
                                   "Formula:\n", 
                                   crayon::green( m_form ) , "\n\n",
                                   "Performance metrics:\n"
                            )
                          )
                          
                          
                          df |>
                            dplyr::select(CV_type, N_folds, N, Random_effect_form, RMSE, Predictive_R2) |>
                            dplyr::rename(`CV type` = CV_type,
                                          `N folds` = N_folds,
                                          # `N pr
                                          `Predictive R2` = Predictive_R2,
                                          `Random-effect structure` = Random_effect_form) |>
                            insight::format_table(use_symbols = T, digits = 0) |>
                            insight::export_table() |>
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
                          
                        #   message(crayon::green(paste0("Saving cross-validation results to: '", output_path, "'\n")) )
                          
                        #   saveRDS(self, output_path)
                        #   assertthat::assert_that(file.exists(output_path),
                        #                           msg = "Writing a copy of the results failed...")
                        #   invisible()
                        # },
                        
                        #' @description
                        #' Function to generate figures summarising the cross-
                        #' validation results.
                        #' 
                        #' @param type (`character()`) Options are: "coefficients",
                        #' "coefficients-summary", "obs-v-pred", "obs-v-error", 
                        #' "pred-v-error", "metrics", "map-sites".
                        #' @param metrics (`character`) Options are the lower-case
                        #' abbreviations of calculated error metrics. - Only needed 
                        #' if type == "metrics".
                        #' @return (`ggplot`).
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
                                              "metrics",
                                              "map-sites"
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
                            
                            orig_coef <- orig_coef |>
                              dplyr::mutate(coef_nm_num = as.numeric(as.factor(Coefficient_name))) |>
                              dplyr::arrange(coef_nm_num)
                            
                            cv_coef <- cv_coef |>
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
                                          y = "Coefficient value",
                                          title = "Variation in fixed-effect coefficients\nover cross-validation folds") +
                              ggplot2::scale_x_continuous(breaks = orig_coef$coef_nm_num, 
                                                          labels = orig_coef$Coefficient_name) +
                              ggplot2::coord_flip() +
                              ggplot2::theme_bw() +
                              basic_thm +
                              ggplot2::facet_wrap(~cv_type) +
                              ggplot2::theme(legend.position = "none")
                            
                          }
                          
                          if (type == "coefficients-summary" ){
                            
                            orig_coef <- self$original_model_obj$fixedef_coefficients |>
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
                                                        
                                                        
                                                        tmp_df <- tmp_df |>
                                                          dplyr::group_by(Coefficient_name, cv_type) |>
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
                            
                            
                            orig_coef <- orig_coef |>
                              dplyr::mutate(coef_nm_num = as.numeric(as.factor(Coefficient_name))) |>
                              dplyr::arrange(coef_nm_num)
                            
                            cv_coef <- cv_coef |>
                              dplyr::ungroup() |>
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
                                            colour = "Cross-\nvalidation\ntype",
                                            title = "Variation in fixed-effect coefficients over cross-validation") +
                              ggplot2::coord_flip() +
                              ggplot2::theme_bw() +
                              basic_thm 
                          }
                          
                          
                          if (type %in% c("obs-v-pred", "obs-v-error", "pred-v-error") ){
                            
                            # Get obs and pred and error - across cv + insample
                            insamp_df <- self$insample_results$insample_pred_df |>
                              dplyr::mutate(cv_type = "In-sample",
                                            fold_id = NA,
                                            Random_effect_form = paste("Random-effects:", Random_effect_form))
                            
                            cv_df <- purrr::map_dfr(.x = self$cv_results,
                                                    .f = function(x){
                                                      # TODO: check nrow(x$cv_pred_df)
                                                      if (is.null(x$cv_pred_df)){
                                                        return(NULL)
                                                      }
                                                      tmp_df <- x$cv_pred_df
                                                      
                                                      cv_type <- x$fold_type
                                                      if (cv_type == "Custom"){
                                                        cv_type <- paste0(cv_type, " (", x$fold_col, ")")
                                                      }
                                                      
                                                      tmp_df$cv_type <- cv_type
                                                      
                                                      tmp_df$fold_id <- as.numeric(as.factor(tmp_df$fold_id))
                                                      
                                                      tmp_df
                                                    })
                            if (nrow(cv_df) > 0){
                              cv_df <- cv_df |>
                              # TODO: only if nrow(cv_df) > 0 do we mutate...
                                dplyr::mutate(Random_effect_form = paste("Random-effects:", Random_effect_form))
                            
                            
                              # Add cv_type that isn't in base_cv_levs to cv_levs
                              add_cv_levs <- unique(cv_df$cv_type)
                              add_cv_levs <- add_cv_levs[!(add_cv_levs %in% base_cv_levs)]
                              cv_levs <- c(base_cv_levs, add_cv_levs)
                              
                              insamp_df$fold_id <- 1
                            } else {
                              message(crayon::yellow("Warning: values from cross-valiation not avaliable."))
                              cv_levs <- base_cv_levs
                            }
                            
                            df <- dplyr::bind_rows(insamp_df,
                                                   cv_df) |>
                              dplyr::mutate(cv_type = factor(cv_type,
                                                             levels = cv_levs)) |>
                              droplevels()
                            
                            
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
                                              y = "Predicted",
                                              title = "Predictions v Observations\nover cross-validation folds") +
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
                                              y = "Error",
                                              title = "Error v Observations\nover cross-validation folds")
                            }
                            if (type == "pred-v-error"){
                              plt <- plt +
                                ggplot2::geom_point(ggplot2::aes(x = prediction,
                                                                 y = error,
                                                                 colour = as.character(fold_id)),
                                                    alpha = 0.25) +
                                ggplot2::labs(x = "Predicted",
                                              y = "Error",
                                              title = "Error v Predictions\nover cross-validation folds")
                            }
                            plt <- plt +
                              ggplot2::theme_bw() +
                              basic_thm +
                              # ggplot2::facet_wrap(~cv_type) +
                              ggplot2::facet_grid(Random_effect_form~cv_type) +
                              
                              ggplot2::theme(legend.position = "none")
                            
                            if (nrow(cv_df) ==  0){
                              plt <- plt +
                                ggplot2::scale_colour_manual(values = c("1" = "#999999"))
                              
                            }
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
                                                      # TODO: check nrow
                                                      if (is.null(x$cv_metrics)){
                                                        return(NULL)
                                                      }
                                                      tmp_df <- x$cv_metrics
                                                      
                                                      cv_type <- x$fold_type
                                                      if (cv_type == "Custom"){
                                                        cv_type <- paste0(cv_type, " (", x$fold_col, ")")
                                                      }
                                                      tmp_df$cv_type <- cv_type
                                                      
                                                      tmp_df
                                                    }) 
                            if (nrow(cv_df) > 0){
                              cv_df <- cv_df |>
                                # TODO: only if nrow(cv_df) > 0 do we filter/mutate...
                                dplyr::filter(Metric_shrt %in% metrics) |>
                                # dplyr::mutate(Random_effect_form = dplyr::case_when(Random_effect_form %in% c("NA", "NULL") ~ paste("Random-effects:", tolower(Random_effect_form)),
                                #                                                     T ~ paste("Random-effects:", Random_effect_form) ))
                                dplyr::mutate(Random_effect_form = paste("Random-effects:\n", Random_effect_form))
                              
                            }
                              
                            
                            
                            # %%%%
                            cv_fold_df <- purrr::map_dfr(.x = self$cv_results,
                                                    .f = function(x){
                                                      # TODO: check nrow
                                                      if (is.null(x$cv_fold_metrics)){
                                                        return(NULL)
                                                      }
                                                      tmp_df <- x$cv_fold_metrics
                                                      
                                                      cv_type <- x$fold_type
                                                      if (cv_type == "Custom"){
                                                        cv_type <- paste0(cv_type, " (", x$fold_col, ")")
                                                      }
                                                      tmp_df$cv_type <- cv_type
                                                      
                                                      tmp_df$fold_id <- as.numeric(tmp_df$fold_id)
                                                      
                                                      tmp_df
                                                    })
                            if (nrow(cv_fold_df) > 0){
                              cv_fold_df <- cv_fold_df |>
                                # TODO: only if nrow(cv_df) > 0 do we filter/mutate...
                                dplyr::filter(Metric_shrt %in% metrics)|>
                                # dplyr::mutate(Random_effect_form = dplyr::case_when(Random_effect_form %in% c("NA", "NULL") ~ paste("Random-effects:", tolower(Random_effect_form)),
                                #                                                     T ~ paste("Random-effects:", Random_effect_form) ))
                                dplyr::mutate(Random_effect_form = paste("Random-effects:\n", Random_effect_form))
                            }
                            
                            # %%%%
                            
                            insamp_df <- self$insample_results$insample_metrics |>
                              dplyr::filter(Metric_shrt %in% metrics)|>
                              # dplyr::mutate(Random_effect_form = dplyr::case_when(Random_effect_form %in% c("NA", "NULL") ~ paste("Random-effects:", tolower(Random_effect_form)),
                              #                                                     T ~ paste("Random-effects:", Random_effect_form) ))
                              dplyr::mutate(Random_effect_form = paste("Random-effects:\n", Random_effect_form))
                            
                            
                            # print(insamp_df)
                            # Add cv_type that isn't in base_cv_levs to cv_levs
                            if (nrow(cv_df) > 0){
                              add_cv_levs <- unique(cv_df$cv_type)
                              add_cv_levs <- add_cv_levs[!(add_cv_levs %in% base_cv_levs)]
                              
                              cv_levs <- base_cv_levs[base_cv_levs %in% unique(cv_df$cv_type)]
                              cv_levs <- c(cv_levs, add_cv_levs)
                              
                              # TODO: check nrow
                              cv_df <- cv_df |>
                                dplyr::mutate(cv_type = factor(cv_type,
                                                               levels = cv_levs),
                                              cv_type_num = as.numeric(cv_type)) |>
                                droplevels()
                              
                              # TODO: check nrow
                              cv_fold_df <- cv_fold_df |>
                                dplyr::mutate(cv_type = factor(cv_type,
                                                               levels = cv_levs),
                                              cv_type_num = as.numeric(cv_type)) |>
                                droplevels()
                            }
                            

                            # TODO: if nrow(cv_df) == 0 - message...
                            if (nrow(cv_df) == 0){
                              message(crayon::yellow(paste0("Warning: metric values from cross-valiation not avaliable.")))
                              
                              plt <- ggplot2::ggplot() +  
                                ggplot2::geom_segment(ggplot2::aes(x = -0.5,
                                                                   xend = +0.5,
                                                                   y = Metric_val,
                                                                   yend = Metric_val),
                                                      data = insamp_df) +
                                ggplot2::geom_rect(ggplot2::aes(xmin = -0.5,
                                                                xmax = +0.5,
                                                                ymin = Metric_lo,
                                                                ymax = Metric_hi),
                                                   alpha = 0.3,
                                                   data = insamp_df) +
                                
                                ggplot2::labs(x = ggplot2::element_blank(),
                                              y = "Metric value",
                                              colour = "Cross-\nvalidation\ntype",
                                              title = "Model performance metrics") +
                                ggplot2::theme_bw() +
                                basic_thm +
                                # ggplot2::facet_wrap(~Metric_fct*Random_effect_form, scales = "free_y",
                                #            labeller = ggplot2::label_parsed) +
                                ggplot2::facet_grid(Metric_fct~Random_effect_form, 
                                                    scales = "free_y",
                                                    labeller = ggplot2::labeller(Metric_fct = ggplot2::label_parsed)
                                                    # labeller = ggplot2::label_parsed
                                ) +
                                ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                               axis.ticks.x = ggplot2::element_blank(),
                                               aspect.ratio = 1)
                              
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
                                
                                cor_df <- cor_df |>
                                  dplyr::filter(Metric %in% insamp_df$Metric) |>
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
                                
                                err_df <- err_df |>
                                  dplyr::filter(Metric %in% insamp_df$Metric) |>
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
                                
                                rse_df <- rse_df |>
                                  dplyr::mutate(Metric_fct = factor(Metric,
                                                                    levels = metr_levs,
                                                                    labels = metr_labs ))
                                
                                plt <- plt +
                                  ggplot2::geom_hline(ggplot2::aes(yintercept = bad), 
                                                      lty = 2, 
                                                      colour = "darkgrey",
                                                      data = rse_df) 
                              }  
                              
                              # stop("No metrics available.",
                              #      call. = FALSE)
                              return(plt)
                            }

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
                              
                              cor_df <- cor_df |>
                                dplyr::filter(Metric %in% cv_df$Metric) |>
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
                              
                              err_df <- err_df |>
                                dplyr::filter(Metric %in% cv_df$Metric) |>
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
                              
                              rse_df <- rse_df |>
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
                                            colour = "Cross-\nvalidation\ntype",
                                            title = "Model performance metrics") +
                              ggplot2::theme_bw() +
                              basic_thm +
                              # ggplot2::facet_wrap(~Metric_fct*Random_effect_form, scales = "free_y",
                              #            labeller = ggplot2::label_parsed) +
                              ggplot2::facet_grid(Metric_fct~Random_effect_form, 
                                                  scales = "free_y",
                                                  labeller = ggplot2::labeller(Metric_fct = ggplot2::label_parsed)
                                                  # labeller = ggplot2::label_parsed
                                                  ) +
                              ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                             axis.ticks.x = ggplot2::element_blank(),
                                             aspect.ratio = 1)
                            
                            
                            
                          }
                          
                          if (type == "map-sites"){
                            # stop("In development")
                            # %%%%%
                            # Just get a coarse-ish world map...
                            tmp <- rnaturalearth::ne_countries(scale = 50, returnclass = 'sf')
                            # NOTE: This differe per model type!!
                            # if (self$model_type == "hib_lme4_model"){
                            #   m_form <- stats::formula(self$original_model_obj$model)
                            # }
                            # if (self$model_type == "hib_inla_model"){
                            #   m_form <- self$original_model_obj$user_formula$full
                            # }
                            # if (self$model_type == "hib_brms_model"){
                            #   m_form <- self$original_model_obj$model$formula$formula
                            # }
                            # 
                            # resp_var <- all.vars(m_form[[2]])
                            # resp_lab <- gsub("_", "\n", resp_var)
                            # 
                            loc_df <- self$original_model_obj$site_locations
                            # %%%
                            
                            plt <- ggplot2::ggplot() +
                              ggplot2::geom_sf(data = tmp,
                                               fill = "lightgrey",
                                               alpha = 0.4,
                                               show.legend = F) +
                              ggplot2::theme_bw() +
                              basic_thm 
                            # %%%
                            
                            cv_df <- purrr::map_dfr(.x = self$cv_results,
                                                     .f = function(x){
                                                       if (is.null(x$cv_pred_df)){
                                                         return(NULL)
                                                       }
                                                       tmp_df <- x$cv_pred_df
                                                       
                                                       cv_type <- x$fold_type
                                                       if (cv_type == "Custom"){
                                                         cv_type <- paste0(cv_type, " (", x$fold_col, ")")
                                                       }
                                                       
                                                       tmp_df$cv_type <- cv_type
                                                       
                                                       tmp_df$fold_id <- as.numeric(as.factor(tmp_df$fold_id))
                                                       
                                                       tmp_df
                                                     })
                            
                            
                            nrow_cv_df <- nrow(cv_df)

                            if (nrow_cv_df > 0){
                              # Add cv_type that isn't in base_cv_levs to cv_levs
                              add_cv_levs <- unique(cv_df$cv_type)
                              add_cv_levs <- add_cv_levs[!(add_cv_levs %in% base_cv_levs)]
                              cv_levs <- c(base_cv_levs, add_cv_levs)
                              
                              cv_df <- cv_df |>
                                dplyr::mutate(cv_type = factor(cv_type,
                                                               levels = cv_levs))
                            }
                            
                            
                            if (nrow_cv_df == 0){
                              # nrow_cv_df <- 0
                              # message(crayon::red("Error: Cross-validation predictions unavailable."))
                              message(crayon::yellow("Warning: map of cross-valiation splits not avaliable."))
                              # return(ggplot2::ggplot() + ggplot2::theme_void())
                              cv_df <- self$insample_results$insample_pred_df |>
                                dplyr::mutate(fold_id = 1,
                                              cv_type = "In-sample",
                                              cv_type = factor(cv_type,
                                                               levels = base_cv_levs))
                              # stop("Cross-validation predictions unavailable.",
                              #      call. = FALSE)
                              # TODO: could use insample df... - all grey...
                            }
                            
                            if ("SSBS" %in% colnames(cv_df)){
                              cv_df <- cv_df |>
                                dplyr::left_join(loc_df, by = "SSBS")
                              
                            }
                            if ("SSBS1" %in% colnames(cv_df)){
                              cv_df <- cv_df |>
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
                            x_rng <- range(cv_df$Longitude, na.rm = T)
                            y_rng <- range(cv_df$Latitude, na.rm = T)
                            
                            
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
                                                               colour = as.factor(as.character(fold_id))),
                                                  alpha = 0.6,
                                                  size = 0.5,
                                                  show.legend = F,
                                                  data = cv_df) +
                              ggplot2::labs(#colour = resp_lab,
                                            x = ggplot2::element_blank(),
                                            y = ggplot2::element_blank(),
                                            title = "Geographic distribution of sites\nand cross-valiation folds") +
                              # ggplot2::scale_color_gradient(low = "firebrick3",
                              #                               high = "dodgerblue3") +
                              ggplot2::xlim(x_rng) +
                              ggplot2::ylim(y_rng) +
                              ggplot2::facet_wrap(~cv_type, ncol = 2)
                            
                            if (nrow_cv_df ==  0){
                              plt <- plt +
                                ggplot2::scale_colour_manual(values = c("1" = "#999999"))
                            }
                            
                          }
                            # %%%%%
                          
                          
                          (plt)
                        }
                      ),
                      private = list()
)
