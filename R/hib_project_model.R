#' hib_project_model: Project a hibiiscus mixed-effects model
#' 
#' @description 
#' This function projects a fitted HibModelXXX using the supplied model object and
#' data to be used in projection.
#' 
#' @details 
#' Given a HibModelXXX, new data and specification, this function projects/predicts 
#' modelled responses for each 'unit' of the new data. 
#' Calculation of the 'Area-of-Applicability' of the model within the new data is 
#' an optional addition.
#' Specification arguments allow for customisation of the prediction type.
#' If specifications are not provided, defaults will be used.
#' Projections will be returned in the same format as the provided data.
#' If data is provided as \code{SpatRaster}, a land-use variable should be included 
#' as a single categorical layer, or as multiple, proportional layers, named to 
#' align with the levels of the categorical variable.
#' Currently, only a single 'categorical' variable can be provided in 
#' 'proportional cover' format in a \code{SpatRaster}.
#' 
#' 
#' @param mod_obj One of; [`HibModelLME4`], [`HibModelINLA`] or [`HibModelBRMS`].
#' @param data A `data.frame`-like object, a \code{SpatRaster}, or a list of 
#' \code{SpatRaster}.
#' Should contain all fixed-effect variables used in the fitted model, and any 
#' required random variables.
#' If `data` is a `data.frame`, `x` and `y` columns should also be included.
#' @param backtransform A `function` that will backtransform the model response
#' to the scale of the original variable, e.g., `exp()` as the inverse of `log()`.
#' Default is `NULL`, i.e., no backtransformation.
#' @param rescale_specs A named `list` (or `NULL`) containing details relevant 
#' to rescaling projections.
#' Specifically, named items should consist of:
#' - `rescale` (`logical`), 
#' - `baseline_covar` (`data.frame`, or NULL), 
#' - `baseline_re_form` (`formula`, `NA`, or `NULL`). \cr
#' 
#' In combination, these arguments determine if predictions will be rescaled or not 
#' (`rescale`), the covariate values used to generate the baseline against which 
#' other predictions are rescaled (`baseline_covar`), and whether/how random 
#' effects should be incorporated in the generation of this baseline (`baseline_re_form`, 
#' see `re_form` in `lme4/inla/brms_specs` argument). \cr
#' If `rescale == TRUE`, but no `baseline_covar` is provided, the intercept of the 
#' model will be used as the baseline value.
#' @param lme4_specs A named `list` (or `NULL`) containing details relevant to 
#' projecting a `HibModelLME4`. \cr
#' Named items should consist of:
#' - `bootstrap` (`logical`), indicating whether bootstrapped prediction intervals 
#' should be estimated. Can be slow for large models.
#' - `n_boot` (`numeric`), how many bootstrap samples should be used to estimate
#' prediction uncertainty.
#' - `calc_se` (`logical`), specify whether or not a standard error based uncertainty 
#' estimate (via [`marginaleffects::predictions`]) is calculated and used to generate prediction 
#' intervals.
#' - `conf_level` (`numeric`) is used to specify the confidence levels of intervals 
#' generated via bootstrapping or the SE-based approach. Values should fall between
#' 0 and 1, typically in the range 0.5 to 0.95.
#' Multiple values can be provided.
#' - `re_form` (`formula`, `NULL`, or `NA`), as `re.form` in [`predict.merMod()`], 
#' used to specify random effects to condition on when predicting.
#' @param inla_specs A named `list` (or `NULL`) containing details relevant to 
#' projecting a `HibModelINLA`. \cr
#' - `re_form` (`formula`, `NULL`, or `NA`), as `re.form` in [`predict.merMod()`], 
#' used to specify random effects to condition on when predicting.
#' - `link_inverse` a `function` that is the inverse of the link used to fit the model.
#' Used to backtranform the global intercept for use in rescaling predictions.
#' @param brms_specs A named `list` (or `NULL`) containing details relevant to 
#' projecting a `HibModelBRMS`.
#' Named items should consist of:
#' - `incl_obs_uncert` (`logical`), to determine whether or not observation-level 
#' uncertainty is included in predictions.
#' Setting this parameter to `TRUE` will produce predictions using
#' [`predict.brmsfit()`] (larger uncertainty), setting to `FALSE` calls 
#' [`fitted.brmsfit()`] internally.
#' - `conf_level` (`numeric`) is used to specify the confidence levels of intervals 
#' generated. Values should fall between 0 and 1, typically in the range 0.5 to 0.95. 
#' Internally converted to quantiles and passed to the `probs` argument of `predict.brmsfit()`/
#' `fitted.brmsfit()`.
#' - `re_form` (`formula`, `NULL`, or `NA`), as `re_formula` in `predict.brmsfit()`, 
#' used to specify random effects to condition on when predicting.
#' @param calc_aoa A `logical` to indicate whether or not the 'Area-of-Applicability' 
#' (AOA) should be estimated - following Meyer and Pebesma (2021).
#' @param cv_obj A [`HibCV`] object, used in AOA estimation.
#' 
#' @return HibProjection
#' 
#' @references 
#' * Meyer, H., & Pebesma, E. (2021). Predicting into unknown space? Estimating the 
#' area of applicability of spatial prediction models. 
#' Methods in Ecology and Evolution, 12, 1620–1633. https://doi.org/10.1111/2041-210X.13650
#' * Meyer H, Milà C, Ludwig M, Linnenbrink J (2023). _CAST: 'caret' Applications 
#' for Spatial-Temporal Models_. R package, <https://CRAN.R-project.org/package=CAST>.
#' 
#' @keywords model, projections
#' 
#' @seealso
#' \code{\link{HibProjection}}
#' \code{\link{hib_calc_aoa}}
#' \code{\link[lme4]{predict.merMod}}
#' \code{\link[brms]{fitted.brmsfit}}
#' \code{\link[brms]{predict.brmsfit}}
#' 
#' @examples 
#' \dontrun{
#'   proj <- hib_mod_project(hib_model, proj_data)
#' } 
#' @export
hib_project_model <- function(mod_obj, 
                              data,
                              backtransform = NULL,
                              rescale_specs = NULL, # list(rescale = T/F, baseline_covar = df, baseline_re_form = formula)
                              lme4_specs = NULL, #list(bootstrap = F, n_boot = NA, calc_se = T, conf_level = 0.95, re_form = NA),
                              inla_specs = NULL, #list(re_form = x, link_inverse = function(y){exp(y)})
                              brms_specs = NULL, #list(incl_obs_uncert = F, conf_level = 0.95, re_form = as lme4...)
                              calc_aoa = FALSE,
                              cv_obj = NULL){

  # TODO: 
  # Better use of functions for model checking/predictions - minor diffs for df v rast

  # Assertions 
  assertthat::assert_that(
    inherits(mod_obj, "HibModelLME4") || inherits(mod_obj, "HibModelINLA") || inherits(mod_obj, "HibModelBRMS"), # obj is a hib model
    is.data.frame(data) || inherits(data, "SpatRaster") || ( is.list(data) & !is.data.frame(data) & all(unlist(lapply(data, function(x){inherits(x, "SpatRaster")}))) ), 
    is.function(backtransform) || is.null(backtransform),
    
    is.list(rescale_specs) || is.null(rescale_specs),

    is.list(lme4_specs) || is.null(lme4_specs),
    is.list(inla_specs) || is.null(inla_specs),
    is.list(brms_specs) || is.null(brms_specs),
    
    is.logical(calc_aoa),
    inherits(cv_obj, "HibCV") || is.null(cv_obj)
  )
  

  # %%%%%
  # Calculation of DI/AOA
  # Need to also incl cv res
  if (calc_aoa){
    assertthat::assert_that(inherits(cv_obj, "HibCV"))  
    message(crayon::green("\nEstimating Area-of-Applicability..."))
    aoa_res <- hib_calc_aoa(mod_obj, cv_obj, data)
  }
  if (!calc_aoa){
    aoa_res <- NULL
  }
  # %%%%%
  
  
  # %%%%%%%%%%
  # Globally, set base_val to NULL, if rescale_pred = F

  # Check/fix rescale_specs
  rescale_specs <- check_proj_specs(rescale_specs, "rescale")
  
  if (rescale_specs$rescale == FALSE){
    baseline_val <- NULL
  }
  # %%%%%%%%%%

  # Check type of model - and process accordingly.
  # %%%%%
  # LME4
  # %%%%%
  if (inherits(mod_obj, "HibModelLME4")){
    message(crayon::green("\nlme4 model detected! Projection in progress..."))
    # TODO: wrap in hib_mod_proj_lme4??
    
    # %%%%%
    # Get info re bootstrapping, re_form, se_calc and conf_level
    # Check/fix lme4 specs
    lme4_specs <- check_proj_specs(lme4_specs, "lme4")
    
    # %%%%%

    message(crayon::green("\nPreparing input data for projection..."))
    # Get some basic info on model covars
    all_covars <- all.vars(stats::formula(stats::delete.response(stats::terms(stats::formula(mod_obj$model)))))
    # But what about random effects... - 
    rand_vars <- names(lme4::ranef(mod_obj$model))
    fix_vars <- all_covars[!all_covars %in% rand_vars]
    
    # TOCHECK - what about random slopes... - how are these presented in model output?

    # then just pass checked lme4_specs to functions - terra_proj...
    # terra_proj_lme4_mod...

    # %%%%%
    # Get baseline value for optional rescaling
    # %%%%%

    # Only needed if rescale_pred == T
    if (rescale_specs$rescale){
      if (is.null(rescale_specs$baseline_covar)){
        # get model intercept
        if ("(Intercept)" %in% mod_obj$fixedef_coefficients$Coefficient_name){
          baseline_val <- mod_obj$fixedef_coefficients$Coefficient_value[mod_obj$fixedef_coefficients$Coefficient_name == "(Intercept)"]
        } else {
          baseline_val <- mod_obj$fixedef_coefficients$Coefficient_value[1]
        }
        # %%%%
        if (inherits(mod_obj$model, "glmerMod")){
          baseline_val <- mod_obj$model@resp$family$linkinv(baseline_val)
        }

        # if backtrans - backrans intercept to get baseline val
        if (!is.null(backtransform)){
          baseline_val <- backtransform(baseline_val)
        }

        baseline_df <- data.frame(Intercept = baseline_val)
      }  

      if (is.data.frame(rescale_specs$baseline_covar)){
        # Check for all fixef and necessary ranefs...
        assertthat::assert_that(all(fix_vars %in% colnames(rescale_specs$baseline_covar)),
                                msg = "Fixed effect variables missing from baseline covariate data.")

        # for categorical vars, ensure levels in projection input are within
        # model fitting data

        check_covar_levels(mod_obj$model, rescale_specs$baseline_covar, fix_vars)

        # need to check any random vars in baseline
        if (inherits(rescale_specs$baseline_re_form, "formula")){
          re_form_vars_baseline <- all.vars(rescale_specs$baseline_re_form)
          # drop any fix effs
          re_form_vars_baseline <- re_form_vars_baseline[!re_form_vars_baseline %in% fix_vars]

          assertthat::assert_that(all(re_form_vars_baseline %in% rand_vars),
                                msg = "Specified re_form structure for baseline prediction contains variables not in the original model.")
          
          check_covar_levels(mod_obj$model, rescale_specs$baseline_covar, re_form_vars_baseline)
        }
        # Add checks if baseline_re_form is NULL.
        if (is.null(rescale_specs$baseline_re_form)){
          check_covar_levels(mod_obj$model, rescale_specs$baseline_covar, rand_vars)
        }

        baseline_val <- stats::predict(mod_obj$model,
                                         rescale_specs$baseline_covar,
                                         type = "response", # lme4_specs$type,
                                         re.form = rescale_specs$baseline_re_form)
        if (!is.null(backtransform)){
          baseline_val <- backtransform(baseline_val)
        }

        baseline_df <- rescale_specs$baseline_covar
        baseline_df$prediction_mean <- baseline_val
      }
    }

    # %%%%%


    # Project for possible data types
    # Dataframe 
    if (is.data.frame(data)){
      # Check all fixef vars are in colnames of data
      assertthat::assert_that(all(fix_vars %in% colnames(data)),
                              msg = "Fixed effect variables missing from input data.")
      
      # for categorical vars, ensure levels in projection input are within 
      # model fitting data
      # for each fix_vars, check if chr or fct
      # if it is, check all data[,var] are in obj$data[,var]
      check_covar_levels(mod_obj$model, data, fix_vars)
      
      # TODO: condense re_form checks
      # need to check any random vars in lme4_specs$re_form
      if (inherits(lme4_specs$re_form, "formula")){
        # print("Checking random variables")
        re_form_vars <- all.vars(lme4_specs$re_form)
        # Drop fixed effects - e.g. from a random slope
        re_form_vars <- re_form_vars[!re_form_vars %in% fix_vars]

        assertthat::assert_that(all(re_form_vars %in% rand_vars), 
                                msg = "Specified re_form structure contains variables not in the original model.")
        # all.vars(lme4_specs$re_form[[2]])
        check_covar_levels(mod_obj$model, data, re_form_vars)
      }
      if (is.null(lme4_specs$re_form)){
        check_covar_levels(mod_obj$model, data, rand_vars)
      }
      
      
      
      # TODO: functionalise the below - for use with terra or df...
      message(crayon::green("\nProjecting model..."))
      pred_df <- as.data.frame(stats::predict(mod_obj$model,
                                       data,
                                       type = "response", # lme4_specs$type,
                                       re.form = lme4_specs$re_form)) 
      colnames(pred_df) <- "prediction_mean" # prediction
      
      # TODO: current approach to backtransforming/rescaling sd doesn't really make sense...

      # %%%
      # Use bootMer to get uncertainty
      if (lme4_specs$bootstrap == T){
        message(crayon::green(paste("\nGenerating bootstrapped prediciton intervals using",
                                lme4_specs$n_boot, "bootstrap samples...")))
        boot_obj <- lme4::bootMer(mod_obj$model, 
                                  function(mod, df = data){
                                    stats::predict(mod, df, 
                                                   type = "response", 
                                                   re.form = lme4_specs$re_form)}, 
                                  nsim = lme4_specs$n_boot,
                                  .progress = "txt")

        # Use conf_level/conf_quants to get prediction_q...
        conf_quants <- conf_level_to_quantiles(lme4_specs$conf_level)

        # add 0.5 to conf_quants, unique values
        conf_quants_ <- unique(c(0.5, conf_quants))
        
        # apply, rename - prediction_boot_q..., prediction_boot_median
        boot_quants <- t(apply(boot_obj$t, 2, function(x) as.numeric(stats::quantile(x, probs = conf_quants_, na.rm = TRUE)))) |> 
          data.frame()

        colnames(boot_quants) <- paste0("prediction_boot_q", conf_quants_*100)
        
        boot_quants <- boot_quants |>
          dplyr::rename(prediction_boot_median = prediction_boot_q50)

        boot_quants$prediction_boot_sd <- apply(boot_obj$t, 2, stats::sd)

        pred_df <- cbind(pred_df, boot_quants)
        # std.err <- apply(boot_test$t, 2, sd)
        # CI.lower1 <- pred - std.err*1.96
        # CI.upper1 <- pred + std.err*1.96
        # %%%
        # Should back trans and rescale here...
      }

      # se-based CIs... and conf level...
      if (lme4_specs$calc_se == T){
        
        message(crayon::green("\nCalculating standard error of model predictions..."))
        # get conf level if specified - only one, or multiple

        # if only 1 conf_level - get single df for se CI
        # prediction_lo_<conf_level>
        se_df <- marginaleffects::predictions(mod_obj$model,
                                              data,
                                              re.form = lme4_specs$re_form,
                                              type = "response",
                                              conf_level = lme4_specs$conf_level[1]) |> 
          data.frame() |> 
          dplyr::select(estimate, std.error, conf.low, conf.high)
        
        colnames(se_df) <- c("estimate",
                             "prediction_se",
                             paste0("prediction_se_", lme4_specs$conf_level[1], "_lo"),
                             paste0("prediction_se_", lme4_specs$conf_level[1], "_hi"))

        # if multiple conf_level - get df once, then use conf level to adjust CI mult and generate se CIs
        if (length(lme4_specs$conf_level) > 1){
          for (i in 2:length(lme4_specs$conf_level)){
            tmp_alpha <- 1 - lme4_specs$conf_level[i]
            tmp_mult <- abs(stats::qnorm(tmp_alpha / 2))
            tmp_lo_nm <- paste0("prediction_se_", lme4_specs$conf_level[i], "_lo")
            tmp_hi_nm <- paste0("prediction_se_", lme4_specs$conf_level[i], "_hi")
            
            se_df[, tmp_lo_nm] <- se_df$estimate - (tmp_mult * se_df$prediction_se)
            se_df[, tmp_hi_nm] <- se_df$estimate + (tmp_mult * se_df$prediction_se)
          }
        }

        se_df <- se_df |>
          dplyr::select(-estimate)

        pred_df <- cbind(pred_df, se_df)
      }

      if (!is.null(backtransform)){
        p_gt_0_orig <- all(pred_df > 0)
        
        pred_df <- backtransform(pred_df)
        
        # - if negative after trans, but none before, print a warning
        p_gt_0_btr <- all(pred_df > 0)
        
        if ((p_gt_0_orig & !p_gt_0_btr)){
          message(crayon::yellow("\nWarning: backtransformation of predictions introduced negative values!"))
        }
      }

      pred_rescaled_df <- NULL
      if (rescale_specs$rescale){
        # divide predictions by baseline val
        pred_rescaled_df <- pred_df / baseline_val #base_val
        colnames(pred_rescaled_df) <- paste("rescaled", colnames(pred_rescaled_df), sep = "_")
      }
      
      data <- dplyr::bind_cols(data, pred_df, pred_rescaled_df)
      
    }
    
    # Single raster object
    if (inherits(data, "SpatRaster")){
      data <- terra_proj(data, mod_obj, backtransform, 
                         rescale_specs$rescale, baseline_val, 
                         fix_vars, rand_vars,
                         proj_specs = lme4_specs)
    }
    # List of raster objects
    if ((is.list(data) & !is.data.frame(data))){
      # & inherits(data[[1]], "SpatRaster")
      # Apply individual terra proj to each item of list...
      data <- purrr::map(data,
                         terra_proj,
                         obj = mod_obj, 
                         backtransform = backtransform, 
                         rescale_pred = rescale_specs$rescale, 
                         base_val = baseline_val, 
                         fix_vars = fix_vars,
                         rand_vars = rand_vars,
                         proj_specs = lme4_specs)
    }
    
  } 
  # %%%%%
  # INLA
  # %%%%%
  else if (inherits(mod_obj, "HibModelINLA")){

    # Package check
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop(
        "Package \"INLA\" must be installed to use this function.",
        call. = FALSE
      )
    }
    # OR
    # rlang::check_installed("INLA", reason = "to use INLA-based model fitting.")
    
    # message(crayon::yellow("Warning: INLA model detected! Projections are not yet supported...\n"))
    message(crayon::green("\nINLA model detected! Projection in progress..."))
    
    # %%%
    # Check specs
    inla_specs <- check_proj_specs(inla_specs, "inla")
    
    # Get some basic info on model covars
    all_covars <- all.vars(stats::formula(stats::delete.response(stats::terms(mod_obj$user_formula$full))))
    # can ignore iid ranefs...
    # but covars are all needed... - can set to NA if missing though
    
    # TODO: put into hib_mod_proj_inla
    
    # Check type of data 
    message(crayon::green("\nPreparing input data for projection..."))
    if (is.data.frame(data)){
      # Don't need to do anything...Just copy data to data_df...
      data_df <- data
    }
    
    #   get from terr_proj_prep...
    # need to know fix_vars for terra_proj style approach...
    fix_vars <- all.vars(stats::formula(stats::delete.response(stats::terms(mod_obj$user_formula$response_fixed))))
    # ignore "Intercept" too
    fix_vars <- fix_vars[!fix_vars == "Intercept"]
    
    # Get random vars
    rand_vars <- all.vars(mod_obj$user_formula$random)
    # drop any vars that are in fix vars... (i.e random slopes..)
    rand_vars <- rand_vars[!rand_vars %in% fix_vars]
    # And re_form_vars
    re_form_vars <- NA
    if (inherits(inla_specs$re_form, "formula")){
      re_form_vars <- all.vars(inla_specs$re_form)
      # drop any vars that are in fix vars... (i.e random slopes..)
      re_form_vars <- re_form_vars[!re_form_vars %in% fix_vars]
      
      assertthat::assert_that(all(re_form_vars %in% rand_vars),
                              msg = "Specified re_form structure contains variables not in the original model.")
    }
    if (is.null(inla_specs$re_form)){
      re_form_vars <- rand_vars
    }


    # Get baseline value for optional rescaling
    # %%%%%
    # If rescale_specs$rescale
    # if rescale_specs$baseline_covar is not a df,
    #  Use intercept
    # Otherwise, check data, as in lme4...

    if (rescale_specs$rescale){
      if (is.null(rescale_specs$baseline_covar)){
        # get model intercept
        int_id <- which(grepl("intercept", rownames(mod_obj$model$summary.fixed), ignore.case = T))
        # What if another covar also contains intercept...
        if (length(int_id) > 1){
          message(crayon::yellow("\nWarning: Multiple 'Intercept' coefficients found, selecting the first..."))
          int_id <- int_id[1]
        }
        # If no explicit intercept, set as 1st value...
        if (length(int_id) == 0){
          int_id <- 1
        }
        baseline_val <- mod_obj$model$summary.fixed$mean[int_id]

        # if backtrans - backrans intercept to get baseline val
        if (!is.null(backtransform)){
          baseline_val <- backtransform(baseline_val)
        }

        baseline_df <- data.frame(Intercept = baseline_val)
      }  

      if (is.data.frame(rescale_specs$baseline_covar)){
        # Check for all fixef and necessary ranefs...
        assertthat::assert_that(all(fix_vars %in% colnames(rescale_specs$baseline_covar)),
                                msg = "Fixed effect variables missing from baseline covariate data.")

        # for categorical vars, ensure levels in projection input are within
        # model fitting data

        check_covar_levels(mod_obj$model, rescale_specs$baseline_covar, fix_vars)

        # need to check any random vars in baseline
        if (inherits(rescale_specs$baseline_re_form, "formula")){
          re_form_vars_baseline <- all.vars(rescale_specs$baseline_re_form)
          # drop any fix effs
          re_form_vars_baseline <- re_form_vars_baseline[!re_form_vars_baseline %in% fix_vars]

          assertthat::assert_that(all(re_form_vars_baseline %in% rand_vars),
                                  msg = "Specified re_form structure for baseline prediction contains variables not in the original model.")
          
          check_covar_levels(mod_obj$model, rescale_specs$baseline_covar, re_form_vars_baseline)
        }

        baseline_df <- rescale_specs$baseline_covar
      }
    }
    # %%%%%

    # If terra, convert rasters to dfs
    if (inherits(data, "SpatRaster")){
      data_df <- terra_to_df_inla(data, mod_obj, fix_vars, rand_vars, re_form_vars)
    }
    # TODO: what about random effects??

    # If list of terra, convert each raster to df, add name of item in list, bind all
    if ((is.list(data) & !is.data.frame(data) &
         all(unlist(lapply(data, function(x){inherits(x, "SpatRaster")}))) ) ){
      data_df <- purrr::map(data,
                            terra_to_df_inla,
                            obj = mod_obj, 
                            fix_vars = fix_vars,
                            rand_vars = rand_vars, 
                            re_form_vars = re_form_vars) |>
        # TODO: better term for id
        dplyr::bind_rows(.id = "Year")
    }
    
    # With data_df...
    # need to check presence of covars, add NA if any missing, add 1 for intercept
    # and lat/long if spatial ranef...
    # obj$inla_fitted_formula
    # all_covars
    # obj$inla_incl_spatial_term
    # obj$inla_stack
    
    
    # %%%%
    
    if (mod_obj$inla_incl_spatial_term){
      # Need s.spde from original model!!
      # or just recreate!
      site.spde <- INLA::inla.spde2.matern(mod_obj$inla_mesh, alpha = 2)
      
      # Need to verify that coords are present in data_df...
      # What if they are missing? - can I use NA
      coord_x_col <- c("longitude",
                       "longitude_midpoint", 
                       "x",
                       "long", 
                       "lon")
      coord_y_col <- c("latitude",
                       "latitude_midpoint", 
                       "y",
                       "lat")
      
      coord_x_col <- colnames(data_df)[tolower(colnames(data_df)) %in% coord_x_col]
      coord_y_col <- colnames(data_df)[tolower(colnames(data_df)) %in% coord_y_col]

      if (length(coord_x_col) == 0){
        message(crayon::yellow("\nWarning: Longitude coordinate column missing! Setting to NA..."))
        coord_x_col <- 'longitude'
        data_df$longitude <- NA
      }
      if (length(coord_y_col) == 0){
        message(crayon::yellow("\nWarning: Latitude coordinate column missing! Setting to NA..."))
        coord_y_col <- 'latitude'
        data_df$latitude <- NA
      }

      assertthat::assert_that(
        length(coord_x_col) == 1,
        length(coord_y_col) == 1
      )

      coords <- as.matrix(data_df[, c(coord_x_col, coord_y_col)])

      s.site <- INLA::inla.spde.make.index('s', n.spde = site.spde$n.spde)
      
      
      if (mod_obj$inla_mesh_settings$mesh_dim == "3d"){
        coords <- INLA::inla.mesh.map(coords,  
                                      projection = "longlat")
      }
      
      pred_A <- INLA::inla.spde.make.A(mod_obj$inla_mesh,
                                       loc = coords)
    } 
    
    # TODO: consider spatial effects in baseline...
    # If coords provided (x, X, Longitude, Longitude_midpoint), 
    #   make an A
    # If not ... predict without space? - how? set pred_A to 0?, have NA coords?

    # do below data prep for baseline_covars
    if (is.data.frame(rescale_specs$baseline_covar)){
      covars_base <- all_covars[all_covars %in% colnames(rescale_specs$baseline_covar)]

      baseline_X <- rescale_specs$baseline_covar |>
        dplyr::select(dplyr::any_of(covars_base))
      if ( ("Intercept" %in% all_covars) &
           (!"Intercept" %in% covars_base)){
        baseline_X$Intercept <- 1
      }
      
      missing_covars <- all_covars[!all_covars %in% colnames(baseline_X)]
      for (covar in missing_covars){
        baseline_X[, covar] <- NA
      }
      missing_covars <- NULL
      
      # for categorical vars, ensure levels in projection input are within 
      # model fitting data (or NA)
      # for each fix_vars, check if chr or fct
      # if it is, check all data[,var] are in obj$data[,var
      check_covar_levels(mod_obj$model, baseline_X, fix_vars)
      
      # also need to check any random vars - difference between fix_vars and all_vars
      if (inherits(rescale_specs$baselin_re_form, "formula")){
          baseline_re_form_vars <- all.vars(rescale_specs$baselin_re_form)
          # drop any vars that are in fix vars... (i.e random slopes..)
          baseline_re_form_vars <- baseline_re_form_vars[!baseline_re_form_vars %in% fix_vars]
      
          # Need to set rand_vars that are non-re_form_vars to NA
          non_reform_rand <- rand_vars[!rand_vars %in% baseline_re_form_vars]
          baseline_X[, non_reform_rand] <- NA
          # check_covar_levels(mod_obj$model, baseline_X, baseline_re_form_vars)
      }
      if (is.na(re_form_vars)){
         # Need to set rand_vars to NA
         baseline_X[, rand_vars] <- NA
      }
      
      if (mod_obj$inla_incl_spatial_term){
        # sort out baseline_A
        # if no coords...

        # stop("Not currently implmented!")

        # %%%%%%
        coord_x_col <- c("longitude",
                         "longitude_midpoint", 
                         "x",
                         "long", 
                         "lon")
        coord_y_col <- c("latitude",
                         "latitude_midpoint", 
                         "y",
                         "lat")
        
        coord_x_col <- colnames(baseline_X)[tolower(colnames(baseline_X)) %in% coord_x_col]
        coord_y_col <- colnames(baseline_X)[tolower(colnames(baseline_X)) %in% coord_y_col]

        if (length(coord_x_col) == 0){
          message(crayon::yellow("\nWarning: Longitude coordinate column missing! Setting to NA..."))
          coord_x_col <- 'longitude'
          baseline_X$longitude <- NA
        }
        if (length(coord_y_col) == 0){
          message(crayon::yellow("\nWarning: Latitude coordinate column missing! Setting to NA..."))
          coord_y_col <- 'latitude'
          baseline_X$latitude <- NA
        }

        assertthat::assert_that(
          length(coord_x_col) == 1,
          length(coord_y_col) == 1
        )

        coords <- as.matrix(baseline_X[, c(coord_x_col, coord_y_col)])

        if (mod_obj$inla_mesh_settings$mesh_dim == "3d"){
          coords <- INLA::inla.mesh.map(coords,  
                                        projection = "longlat")
        }
        
        baseline_A <- INLA::inla.spde.make.A(mod_obj$inla_mesh,
                                       loc = coords)
        # %%%%%%

        baseline_stk <- INLA::inla.stack(tag = "baseline",
                                         data = list(y = rep(NA, nrow(baseline_X))),
                                         A = list(baseline_A, 1),
                                         effects = list(s = s.site,
                                                        X = baseline_X))
      } else {
        baseline_stk <- INLA::inla.stack(tag = "baseline",
                                     data = list(y = rep(NA, nrow(baseline_X))),
                                     A = list(1),
                                     effects = list(X = baseline_X))
      }

    } 

    # # Set up X data for modelling
    # pred_X <- All cols in data-df, apart form lat/long
    # if Intercept missing, add
    # if any other cols missing, e.g. SS/SSB, add as NA
    covars_sub <- all_covars[all_covars %in% colnames(data_df)]
    

    pred_X <- data_df |>
      dplyr::select(dplyr::any_of(covars_sub))
    
    if ( ("Intercept" %in% all_covars) &
         (!"Intercept" %in% covars_sub)){
      pred_X$Intercept <- 1
    }
    
    missing_covars <- all_covars[!all_covars %in% colnames(pred_X)]
    for (covar in missing_covars){
      pred_X[, covar] <- NA
    }
    
    # for categorical vars, ensure levels in projection input are within 
    # model fitting data (or NA)
    # for each fix_vars, check if chr or fct
    # if it is, check all data[,var] are in obj$data[,var
    check_covar_levels(mod_obj$model, pred_X, fix_vars)
    
    # also need to check any random vars - difference between fix_vars and all_vars
    if (inherits(inla_specs$re_form, "formula")){
        # Need to set rand_vars that are non-re_form_vars to NA
        non_reform_rand <- rand_vars[!rand_vars %in% re_form_vars]
        pred_X[, non_reform_rand] <- NA
        check_covar_levels(mod_obj$model, pred_X, re_form_vars)
    }
    if (is.null(inla_specs$re_form)){
        check_covar_levels(mod_obj$model, pred_X, re_form_vars)
    }
    if (is.na(re_form_vars)){
       # Need to set rand_vars to NA
       pred_X[, rand_vars] <- NA
    }

    # Create prediction stack
    if (mod_obj$inla_incl_spatial_term){
      pred_stk <- INLA::inla.stack(tag = "pred",
                                   data = list(y = rep(NA, nrow(pred_X))),
                                   A = list(pred_A, 1),
                                   effects = list(s = s.site,
                                                  X = pred_X))
    } else {
      pred_stk <- INLA::inla.stack(tag = "pred",
                                   data = list(y = rep(NA, nrow(pred_X))),
                                   A = list(1),
                                   effects = list(X = pred_X))
    }
    # Get original model fitting stack..
    # And combine with pred for refitting

    # if baseline_covar, add in too
    if (is.data.frame(rescale_specs$baseline_covar)){
      comb_stk <- INLA::inla.stack(mod_obj$inla_stack,
                                   pred_stk,
                                   baseline_stk)
      } else {
      comb_stk <- INLA::inla.stack(mod_obj$inla_stack,
                                   pred_stk)
    }

    
    message(crayon::green("\nRe-fitting INLA model to make projection (this may take a while)..."))
    # re-fit/ Project
    # incl while to ensure fixef are fitted properly!
    # %%%%%
    if (length(fix_vars) > 0){
      check_flag <- 0
      while(check_flag == 0){
        m0 <- INLA::inla(formula = mod_obj$inla_fitted_formula,
                   family = mod_obj$model$.args$family,
                   data = INLA::inla.stack.data(comb_stk),
                   control.fixed = list(expand.factor.strategy='inla'),
                   # control.inla = list(strategy = "", int.strategy = ""), 
                   control.compute = list(dic = TRUE, waic = TRUE),
                   control.predictor = list(A = INLA::inla.stack.A(comb_stk),
                                            link = 1),
                   verbose = F)

        check_flag <- nrow(m0$summary.fixed)
      }
    } else {
      m0 <- INLA::inla(formula = mod_obj$inla_fitted_formula,
                 family = mod_obj$model$.args$family,
                 data = INLA::inla.stack.data(comb_stk),
                 control.fixed = list(expand.factor.strategy='inla'),
                 # control.inla = list(strategy = "", int.strategy = ""), 
                 control.compute = list(dic = TRUE, waic = TRUE),
                 control.predictor = list(A = INLA::inla.stack.A(comb_stk),
                                          link = 1),
                 verbose = F)
    }
    # %%%%%

    # Get pred_stack index
    pred_idx <- INLA::inla.stack.index(comb_stk, tag = "pred")$data
    
    pred_df <- m0$summary.fitted.values[pred_idx, ]
    rownames(pred_df) <- NULL
    
    colnames(pred_df) <- paste("prediction", colnames(pred_df), sep = "_")
    
    # Harmonise prediction colnames
    
    # prediction_mean prediction_sd prediction_q2.5
    # prediction_0.5quant prediction_0.975quant prediction_mode
    
    # prediction_mean prediction_sd prediction_0.025quant
    # prediction_median prediction_q97.5 prediction_mode
    pred_df <- pred_df |>
      dplyr::rename(prediction_q2.5   = prediction_0.025quant,
                    prediction_median = prediction_0.5quant,
                    prediction_q97.5  = prediction_0.975quant)
    
    if (!is.null(backtransform)){
      p_gt_0_orig <- all(pred_df > 0)
      pred_df <- backtransform(pred_df)
      p_gt_0_btr <- all(pred_df > 0)
      if ((p_gt_0_orig & !p_gt_0_btr)){
        message(crayon::yellow("\nWarning: backtransformation of predictions introduced negative values!"))
      }
    }
    
    pred_df <- dplyr::bind_cols(data_df, pred_df)
    
    # Need to aggr proportional cover before rescaling...
    # But only if raster??
    if ((inherits(data, "SpatRaster")) |
        (is.list(data) & !is.data.frame(data) & 
          all(unlist(lapply(data, function(x){inherits(x, "SpatRaster")}))) ) ){
      if ("prop_cover" %in% colnames(data_df)){
        # How to better deal with "Year" - could be sc_yr
        if ("Year" %in% colnames(data_df)){
          pred_df <- pred_df |>
            dplyr::group_by(x, y, Year) |>
            dplyr::summarise(dplyr::across(dplyr::contains("prediction"), ~ stats::weighted.mean(.x, w = prop_cover)))
        }
        if (!"Year" %in% colnames(data_df)){
          pred_df <- pred_df |>
            dplyr::group_by(x, y)  |>
            dplyr::summarise(dplyr::across(dplyr::contains("prediction"), ~ stats::weighted.mean(.x, w = prop_cover)))
          # weighted av of all predictions cols  
        }
      }
    }
    
    # Optional rescaling
    pred_rescaled_df <- NULL
    
    if (rescale_specs$rescale){
      if (is.data.frame(rescale_specs$baseline_covar)){
        # Get baseline
        baseline_idx <- INLA::inla.stack.index(comb_stk, tag = "baseline")$data

        baseline_val <- m0$summary.fitted.values[baseline_idx, ]
        rownames(baseline_val) <- NULL
    
        colnames(baseline_val) <- paste("prediction", colnames(baseline_val), sep = "_")
    
        # Harmonise prediction colnames
        
        # prediction_mean prediction_sd prediction_q2.5
        # prediction_0.5quant prediction_0.975quant prediction_mode
        
        # prediction_mean prediction_sd prediction_0.025quant
        # prediction_median prediction_q97.5 prediction_mode
        baseline_val <- baseline_val |>
          dplyr::rename(prediction_q2.5   = prediction_0.025quant,
                        prediction_median = prediction_0.5quant,
                        prediction_q97.5  = prediction_0.975quant)
        # Backtrans
        if (!is.null(backtransform)){
          baseline_val <- backtransform(baseline_val)
        }
        # Add to df
        baseline_df <- cbind(baseline_df, baseline_val)
        baseline_val <- baseline_df$prediction_mean
      }

      # Otherwise, already have baseline_val (intercept)

      # divide predictions by baseline val
      pred_cols <- colnames(pred_df)[grepl("prediction_", colnames(pred_df))]
      pred_rescaled_df <- pred_df[,pred_cols] / baseline_val
      colnames(pred_rescaled_df) <- paste("rescaled", colnames(pred_rescaled_df), sep = "_")
    }
    
    pred_df <- dplyr::bind_cols(pred_df, pred_rescaled_df)
    
    # make back into data... - if df
    if (is.data.frame(data)){
      data <- pred_df
    }
    pred_cols <- colnames(pred_df)[grepl("prediction_", colnames(pred_df))]
    
    # if input was terra/list of terra, convert pred df back to terra format
    if (inherits(data, "SpatRaster")){
      data_pred <- terra::rast(pred_df[,c("x", "y", pred_cols)],
                               crs = terra::crs(data))
      
      data <- c(data, data_pred)
    }
    if ((is.list(data) & !is.data.frame(data) &
         all(unlist(lapply(data, function(x){inherits(x, "SpatRaster")}))) )){
      data <- purrr::map(.x = unique(pred_df$Year),
                         function(.x, df = pred_df, 
                                  data_orig = data, 
                                  col_nms = pred_cols){
                           tmp_df <- df |>
                             dplyr::filter(Year == .x) |>
                             dplyr::select(-Year)
                           
                           tmp_df <- tmp_df[,c("x", "y", col_nms)]
                           
                           tmp_r <- terra::rast(tmp_df,
                                                crs = terra::crs(data_orig[[1]]))
                           
                           out <- c(data_orig[[.x]], tmp_r)
                         })
      
      names(data) <- as.character(unique(pred_df$Year))
    }
    
    
    # %%%%%
    # BRMS
    # %%%%%
  } else if (inherits(mod_obj, "HibModelBRMS")){
    message(crayon::green("\nbrms model detected! Projection in progress..."))
    
    # %%%%
    # Check specs
    brms_specs <- check_proj_specs(brms_specs, "brms")
    
    
    message(crayon::green("\nPreparing input data for projection..."))
    # Get some basic info on model covars
    all_covars <- all.vars(stats::formula(stats::delete.response(stats::terms(mod_obj$model$formula$formula))))
    rand_vars <- names(brms::ranef(mod_obj$model))
    fix_vars <- all_covars[!all_covars %in% rand_vars]

    # %%%%%
    # Only needed if rescale_pred == T
    if (rescale_specs$rescale){
      if (is.null(rescale_specs$baseline_covar)){
        # get model intercept
        if ("Intercept" %in% mod_obj$fixedef_coefficients$Coefficient_name){
          baseline_val <- mod_obj$fixedef_coefficients$Coefficient_value[mod_obj$fixedef_coefficients$Coefficient_name == "Intercept"]
        } else {
          baseline_val <- mod_obj$fixedef_coefficients$Coefficient_value[1]
        }
        # %%%%
        
        baseline_val <- mod_obj$model$family$linkinv(base_val)
        

        # if backtrans - backrans intercept to get baseline val
        if (!is.null(backtransform)){
          baseline_val <- backtransform(baseline_val)
        }

        baseline_df <- data.frame(Intercept = baseline_val)
      }  

      if (is.data.frame(rescale_specs$baseline_covar)){
        # Check for all fixef and necessary ranefs...
        assertthat::assert_that(all(fix_vars %in% colnames(rescale_specs$baseline_covar)),
                                msg = "Fixed effect variables missing from baseline covariate data.")

        # for categorical vars, ensure levels in projection input are within
        # model fitting data
        check_covar_levels(mod_obj$model, rescale_specs$baseline_covar, fix_vars)

        # need to check any random vars in baseline
        if (inherits(rescale_specs$baseline_re_form, "formula")){
          re_form_vars_baseline <- all.vars(rescale_specs$baseline_re_form)
          # drop any fix effs
          re_form_vars_baseline <- re_form_vars_baseline[!re_form_vars_baseline %in% fix_vars]

          assertthat::assert_that(all(re_form_vars_baseline %in% rand_vars),
                                msg = "Specified re_form structure for baseline prediction contains variables not in the original model.")
          
          check_covar_levels(mod_obj$model, rescale_specs$baseline_covar, re_form_vars_baseline)
        }
        if (is.null(rescale_specs$baseline_re_form)){
          check_covar_levels(mod_obj$model, rescale_specs$baseline_covar, rand_vars)
        }

        baseline_val <- stats::fitted(mod_obj$model,
                                         rescale_specs$baseline_covar,
                                         scale = "response",
                                         re_form = rescale_specs$baseline_re_form) |>
          as.data.frame()|>
          dplyr::rename(prediction_mean  = Estimate,
                        prediction_sd    = Est.Error,
                        prediction_q2.5  = Q2.5,
                        prediction_q97.5 = Q97.5)

        if (!is.null(backtransform)){
          baseline_val <- backtransform(baseline_val)
        }

        baseline_df <- rescale_specs$baseline_covar
        baseline_df <- cbind(baseline_df, baseline_val)

        baseline_val <- baseline_val$prediction_mean
      }
    }

    # %%%%%
    
    # Project for possible data types
    # Dataframe 
    if (is.data.frame(data)){
      # Check all fixef vars are in colnames of data
      assertthat::assert_that(all(fix_vars %in% colnames(data)),
                              msg = "Fixed effect variables missing from input data.")
      
      
      # for categorical vars, ensure levels in projection input are within 
      # model fitting data
      # for each fix_vars, check if chr or fct
      # if it is, check all data[,var] are in obj$data[,var
      check_covar_levels(mod_obj$model, data, fix_vars)
      
      # also need to check any random vars in brms_specs$re_form
      # Check they are present, and check they have the correct levels...
      if (inherits(brms_specs$re_form, "formula")){
          re_form_vars <- all.vars(brms_specs$re_form)
          re_form_vars <- re_form_vars[!re_form_vars %in% fix_vars]
          
          assertthat::assert_that(all(re_form_vars %in% rand_vars), 
                                  msg = "Specified re_form structure contains variables not in the original model.")
          check_covar_levels(mod_obj$model, data, re_form_vars)
      }
      if (is.null(brms_specs$re_form)){
        check_covar_levels(mod_obj$model, data, rand_vars)
      }

      # Convert conf levels to quantiles
      conf_quants <- conf_level_to_quantiles(brms_specs$conf_level)

      # %%%%%
      message(crayon::green("\nProjecting model..."))
      if (brms_specs$incl_obs_uncert == TRUE){
      
        mean_df <- stats::predict(mod_obj$model,
                           data,
                           re_form = brms_specs$re_form,
                           probs = conf_quants) |>
          as.data.frame()
        med_df <- stats::predict(mod_obj$model,
                          data,
                          re_form = brms_specs$re_form,
                          robust = TRUE) |>
          as.data.frame()
      }
      if (brms_specs$incl_obs_uncert == FALSE){
        # brms::fitted.brmsfit
        mean_df <- stats::fitted(mod_obj$model,
                          data,
                          scale = "response", # brms_specs$scale,
                          re_form = brms_specs$re_form,
                           probs = conf_quants) |>
          as.data.frame()
          # brms::fitted.brmsfit
        med_df <- stats::fitted(mod_obj$model,
                         data,
                         scale = "response",# brms_specs$scale,
                         re_form = brms_specs$re_form,
                         robust = TRUE) |>
          as.data.frame()
      }

      mean_df <- mean_df |>
        dplyr::rename(prediction_mean  = Estimate,
                      prediction_sd    = Est.Error) |>
        dplyr::rename_with(~gsub("Q", "prediction_q", .x, fixed = T))#,

      med_df <- med_df |>
        dplyr::rename(prediction_median= Estimate,
                      prediction_mad   = Est.Error)
      
      pred_df <- mean_df |>
        dplyr::bind_cols(med_df |>
                           dplyr::select(prediction_median,
                                         prediction_mad))
      # %%%%%

      if (!is.null(backtransform)){
        p_gt_0_orig <- all(pred_df > 0)

        pred_df <- backtransform(pred_df)
        
        p_gt_0_btr <- all(pred_df > 0)

        if ((p_gt_0_orig & !p_gt_0_btr)){
          message(crayon::yellow("\nWarning: backtransformation of predictions introduced negative values!"))
        }
      }
      # %%%
      # Optional rescaling
      pred_rescaled_df <- NULL

      if (rescale_specs$rescale){
        # divide predictions by baseline val
        pred_cols <- colnames(pred_df)[grepl("prediction_", colnames(pred_df))]
        pred_rescaled_df <- pred_df[,pred_cols] / baseline_val
        colnames(pred_rescaled_df) <- paste("rescaled", colnames(pred_rescaled_df), sep = "_")
      }
      # pred_rescaled_df <- pred_rescaled_df |>
      #   dplyr::select(-rescaled_prediction_error)
      
      data <- dplyr::bind_cols(data, pred_df, pred_rescaled_df)
    }
    
    # Single raster object
    if (inherits(data, "SpatRaster")){
      data <- terra_proj(data, mod_obj, backtransform, rescale_specs$rescale, baseline_val, fix_vars,
                          rand_vars,
                         proj_specs = brms_specs)
    }
    # List of raster objects
    if (( is.list(data) & !is.data.frame(data) & 
          all(unlist(lapply(data, function(x){inherits(x, "SpatRaster")}))) ) ){
      # Apply individual terra proj to each item of list...
      data <- purrr::map(data,
                         terra_proj,
                         obj = mod_obj,
                         backtransform = backtransform,
                         rescale_pred = rescale_specs$rescale, # rescale_pred,
                         base_val = baseline_val, # base_val,
                         fix_vars = fix_vars,
                         rand_vars = rand_vars,
                         proj_specs = brms_specs)
    }
    
  }
  
  # create projection object - need model obj, projection output, and projection_specs (backrans, rescale)

  if (rescale_specs$rescale & 
    is.data.frame(rescale_specs$baseline_covar)){
    # update baseline_covar to baseline_df
    rescale_specs$baseline_df <- baseline_df
    # don't need baseline_covar??
  } else {
    rescale_specs$baseline_df <- NA
  }
  
  # TODO: flag aoa calc

  proj <- HibProjection$new(obj = mod_obj,
                                 proj = data,
                                 backtransform = backtransform,
                                 rescale_specs = rescale_specs,
                                 lme4_specs = lme4_specs,
                                 inla_specs = inla_specs,
                                 brms_specs = brms_specs,
                                 aoa_res    = aoa_res)
  return(proj)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal helper functions for projection
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' terra_proj: A helper function to wrap the preparation of SpatRaster data, and 
#' then project a hib, mixed-effects model - for lme4 and brms models only.
#' 
#' @description 
#' This function projects a fitted hib model using the supplied model object and
#' `SpatRaster` data to be used in projection.
#' 
#' @param data A `SpatRaster` that can be used in model projection/prediction.
#' @param obj One of; `HibModelLME4` or `HibModelBRMS`.
#' @param backtransform A `function` that will backtransform the model response
#' to the scale of the original variable, e.g., `exp()` as the inverse of `log()`.
#' @param rescale_pred A `logical` value - should the predictions be rescaled 
#' relative to a baseline value.
#' @param base_val A baseline value used in rescaling predictions.
#' @param fix_vars A `character` vector of model fixed effects.
#' @param rand_vars A `character` vector of model random effects.
#' @param proj_specs A `list` containing details relevant to projecting the model,
#' either lme4_specs, or brms_specs.
#' 
#' @return SpatRaster
#' 
#' @keywords model, projections
#' 
#' @seealso 
#' \code{\link{hib_project_model}}
#' 
#' @examples 
#' \dontrun{
#'   tmp_proj <- terra_proj(data, obj, backtransform, rescale_pred, base_val, fix_vars, rand_vars, proj_specs)
#' } 
#' @noRd
terra_proj <- function(data, 
                       obj, 
                       backtransform, 
                       rescale_pred, 
                       base_val, 
                       fix_vars,
                       rand_vars,
                       proj_specs = NULL
                       ){
  
  # Assertions
  assertthat::assert_that(
    inherits(obj, "HibModelLME4") || inherits(obj, "HibModelBRMS"),
    inherits(data, "SpatRaster"), 
    is.function(backtransform) || is.null(backtransform),
    is.logical(rescale_pred),
    is.numeric(base_val) || is.null(base_val),
    is.character(fix_vars),
    is.character(rand_vars) || is.na(c(rand_vars)) || is.null(rand_vars),
    is.list(proj_specs) || is.null(proj_specs)
  )
 
  # Extract important re_form_vars
  # Condense into single set of calls to specs
  re_form_vars <- NA
  if (inherits(proj_specs$re_form, "formula")){
    re_form_vars <- all.vars(proj_specs$re_form)
    re_form_vars <- re_form_vars[!re_form_vars %in% fix_vars]
    assertthat::assert_that(all(re_form_vars %in% rand_vars), 
                            msg = "Specified re_form structure contains variables not in the original model.")
  }
  if (is.null(proj_specs$re_form)){
    re_form_vars <- rand_vars
  }
  # For categorical LUI rasters
  if (all(fix_vars %in% names(data))){
    # Covar checks
    check_covar_levels(obj$model, data, fix_vars)

    if ( (inherits(proj_specs$re_form, "formula") | 
          is.null(proj_specs$re_form)) ){
      check_covar_levels(obj$model, data, re_form_vars)
    }
    
    
    if (inherits(obj, "HibModelLME4")){
      data_pred <- terra_proj_lme4_pred(obj$model, data, backtransform, proj_specs)
    }

    if (inherits(obj, "HibModelBRMS")){
      data_pred <- terra_proj_brms_pred(obj$model, data, backtransform, proj_specs)
    }
    
    # %%%%
  } else {
    # Categorical variable present as proportional cover ...
    # Any layer with a name in fix_vars is good to use 
    covars_in_rast <- names(data)[names(data) %in% fix_vars]
    # incl ranefs here -- either from re_form, or model
    if (!all(is.na(re_form_vars))){
      ranefs_in_rast <- names(data)[names(data) %in% re_form_vars]
      covars_in_rast <- c(covars_in_rast, ranefs_in_rast)
    }

    tmp_pred_rast <- NULL
    n_pred_rast_lyrs <- 0
    if (length(covars_in_rast) > 0 ){
      tmp_pred_rast <- data[[covars_in_rast]]  
      n_pred_rast_lyrs <- dim(tmp_pred_rast)[3]
    }
    
    # any layer with a name not in fix_covars needs to be aligned to a var
    non_covars_in_rast <- names(data)[!names(data) %in% fix_vars]
    # drop ranefs here -- either from re_form, or model ??
    # Easy way to drop ranefs_in_rast 
    non_covars_in_rast <- non_covars_in_rast[!non_covars_in_rast %in% covars_in_rast]
    # But, also want to drop ranefs that aren't in rast from non_covars_in_rast
    non_covars_in_rast <- non_covars_in_rast[!non_covars_in_rast %in% rand_vars]

    tmp_rast <- NULL
    if (length(non_covars_in_rast) > 0 ){
      tmp_rast <- data[[non_covars_in_rast]]
    }
    
    # Assume that layer names match levels of covariates
    # e.g "LUIPrimary_minimal use" or just "Primary_minimal use"
    covars_not_in_rast <- fix_vars[!fix_vars %in% names(data)]
    
    if (length(non_covars_in_rast) > 0 ){
      # for each covar, get unique levels (if factor in mod_df) - what about chr..
      # Use model-fitting data from obj$model...
      # lmer: @frame
      # inla: $.args$data
      # brms: @frame
      
      # lme4 fit data,
      if (inherits(obj$model, "merMod")){
        fit_dat <- obj$model@frame
      }
      # brms fit data
      if (inherits(obj$model, "brmsfit")){
        fit_dat <- obj$model$data
      }
      # inla fit data - not aplicable
      # if (class(obj$model) == "inla"){
      # if (inherits(obj$model, "inla")){
      #   fit_dat <- obj$model$.args$data
      # }
      # 
      covars_not_in_rast <- purrr::map_dfr(covars_not_in_rast,
                                           function(x){
                                             if (is.character(fit_dat[[x]])){
                                               tmp_out <- data.frame(covar_nm = x,
                                                                     covar_levs = unique(fit_dat[[x]]))
                                               # Could add covar_lev_idx here... - save doing later...
                                             } else if (is.factor(fit_dat[[x]])){
                                               tmp_out <- data.frame(covar_nm = x,
                                                                     covar_levs = levels(fit_dat[[x]]))
                                               # Could add covar_lev_idx here... - save doing later...
                                             } else {
                                               tmp_out <- data.frame(covar_nm = x,
                                                                     covar_levs = NA)
                                             }
                                             tmp_out
                                           })
      
      covars_not_in_rast <- covars_not_in_rast |>
        # Combine covar nm with levels
        dplyr::mutate(covar_nm_levs = paste0(covar_nm, covar_levs), 
                      # check for presence in tmp_rast
                      covar_levs_in_rast = covar_levs %in% names(tmp_rast), 
                      covar_nm_levs_in_rast = covar_nm_levs %in% names(tmp_rast),
                      # Add level name in raster
                      covar_levs_rast = dplyr::case_when(covar_levs_in_rast == T ~ covar_levs,
                                                         covar_nm_levs_in_rast == T ~ covar_nm_levs,
                                                         T ~ NA))
      
      # Check for binary layers... 
      # TODO: what about NA.../nan... (need to account for NA values...)
      is_binary_layer <- purrr::map_dfr(names(tmp_rast),
                                        function(nm){
                                          data.frame(covar_levs_rast = nm,
                                                     is_binary_layer = all(unique(terra::values(tmp_rast[[nm]])) %in% c(0,1, NA, NaN)) # NA, NaN
                                          )
                                        })
      
      covars_not_in_rast <- covars_not_in_rast |>
        dplyr::left_join(is_binary_layer, by = "covar_levs_rast") |>
        dplyr::group_by(covar_nm) |>
        dplyr::mutate(is_binary_covar = all(is_binary_layer, na.rm = T))
      
      
      # Split binary and non binary layers
      bin_covars <- covars_not_in_rast |>
        dplyr::filter(is_binary_covar == T)
      
      non_bin_covars <- covars_not_in_rast |>
        dplyr::filter(is_binary_covar == F)
      
      # TODO: maybe this can be dropped???
      # If we don't assume any within cell spatial structure, just sum/multiply
      # proportional effects of each covariate level...
      # BUT, this approach doesn't work for terra::predict
      # would have to edit.../predict using only prop vars,
      # then do multiplication to aggregate
      # predict for non-prop vars - sum
      if (length(unique(non_bin_covars$covar_nm)) > 1){
        message(crayon::yellow("\nWarning: SpatRaster with multiple non-binary layers detected! Projections are not yet supported..."))
        stop("In development")
      }
      
      
      # %%%%%%%%%%%%%%%%%%%%
      # If there is a non_covar in rast that doesn't have a matching covar/lev
      # throw an error.
      # %%%%%%%%%%%%%%%%%%%%
      
      # ID layers that are not within covariates required
      problem_layers <- non_covars_in_rast[!non_covars_in_rast %in% covars_not_in_rast$covar_levs_rast]
      
      # Are there other situations where a problem layer could arise - or actualy not be a problem??
      # Note: extra layers in input that are not relevant to the model - eg a non-existent
      # covar level (Not Pasture) is simply ignored!
      if (length(problem_layers) > 0){
        # Maybe better to throw a warning/error??
        message(crayon::red(paste0("\nError: covariate layer/s in the projection ",
                                   "data do not align with the model fitting data!\n\n",
                                   "Please fix this, see details below:\n")
                            )
                )      
  
        message(paste("Problematic layers: ",
                      paste(problem_layers, collapse = " ; "),
                      "\n") 
                )
        
        stop("Unknown covariate layer, see above.", call. = F)
        
        # TODO: Could I just drop these layers? - e.g, area layer in globiom data
      }
      
      
      # %%%%%%%%%%%%%%%%%%%%%
      
      # For binary layers - factor is layer name, minus covar nm if name is covar_nm_lev
      # Collapse to single stack with layer names as covar_nm
      
      # only if binary layers
      tmp_bin_rast <- NULL
      n_bin_covar_layers <- 0
      if (nrow(bin_covars) > 0){
        tmp_bin_rast <- purrr::map(unique(bin_covars$covar_nm),
                                   terra_covar_prep,
                                   covar_df = bin_covars,
                                   data_rast = tmp_rast,
                                   mod_data = fit_dat,
                                   input_type = "binary") |>
          unlist() |>
          c() |>
          terra::rast()
        
        n_bin_covar_layers <- dim(tmp_bin_rast)[3]
      }
      
      # For non-binary layers - factor is layer name, minus covar nm if name is covar_nm_lev
      # BUT, also need to record values for later aggregation
      # currently assume only one non_bin_covar...
      
      tmp_prop_rast <- NULL
      n_prop_covar_lyrs <- 0
      if (nrow(non_bin_covars) > 0){
        # No need for purrr - only gets here if 1 non-binary covar...
        tmp_prop_rast <- purrr::map(unique(non_bin_covars$covar_nm),
                                    terra_covar_prep,
                                    covar_df = non_bin_covars,
                                    data_rast = tmp_rast,
                                    mod_data = fit_dat,
                                    input_type = "non-binary") |>
          unlist() |>
          c() |>
          terra::rast()
        # Want to end up with N raster stacks - 1 per prop level
        n_prop_covar_lyrs <- dim(tmp_prop_rast)[3]
      }
      # Check tmp_prop_rast is a terra SpatRaster
      # Could optionally exit here if prepping for inla...
      
    }
    
    # check N layers of each
    # if n_prop_covar_lyrs > 1  - list of rasters - each named based on covar level
    if (n_prop_covar_lyrs > 1){
      
      non_bin_sub <- non_bin_covars |>
        dplyr::filter(!is.na(covar_levs_rast))
      
      non_bin_r <- tmp_rast[[non_bin_sub$covar_levs_rast]]
      
      tmp_pred_rast_ls <- vector("list", n_prop_covar_lyrs)
      names(tmp_pred_rast_ls) <- rep("", n_prop_covar_lyrs)

      # Optionally include se message 
      # Then per iter, print lyr i...

      for (i in 1:n_prop_covar_lyrs){
        # TODO: need to be careful combining these, if first item is null,
        # output becomes a list of 1...
        
        tmp_ <- c(tmp_pred_rast, tmp_bin_rast, tmp_prop_rast[[i]])
        # If tmp is a list, get [[1]]
        if (is.list(tmp_)){
          tmp_ <- tmp_[[1]]
        }

        assertthat::assert_that(
          inherits(tmp_, "SpatRaster")
        )
        # If is a SpatRaser, use as is.

        tmp_pred_rast_ls[[i]] <- tmp_ 
        
        # Check all covars in new SpatRaster
        assertthat::assert_that(all(fix_vars %in% names(tmp_pred_rast_ls[[i]])))
        
        check_covar_levels(obj$model, tmp_pred_rast_ls[[i]], fix_vars)
        
        if ((inherits(proj_specs$re_form, "formula") |
             is.null(proj_specs$re_form)) ){
          check_covar_levels(obj$model, tmp_pred_rast_ls[[i]], re_form_vars)
        }
        
        if (inherits(obj, "HibModelLME4")){
          # TODO: silent proj
          tmp_pred_rast_ls[[i]] <- terra_proj_lme4_pred(obj$model, 
                                                        tmp_pred_rast_ls[[i]], 
                                                        backtransform,
                                                        proj_specs)
        }
        if (inherits(obj, "HibModelBRMS")){
          # TODO: silent proj
          tmp_pred_rast_ls[[i]] <- terra_proj_brms_pred(obj$model, 
                                                        tmp_pred_rast_ls[[i]], 
                                                        backtransform,
                                                        proj_specs)
        }
        
        
        # Get proportional weighting
        lyr_nm <- terra::cats(tmp_prop_rast)[[i]] |> 
          dplyr::filter(ID == as.vector(unname(unique(terra::values(tmp_prop_rast[[i]])))) ) |> 
          dplyr::pull(category)
        
        lyr_orig_nm <- non_bin_sub |>
          dplyr::filter(covar_levs == lyr_nm) |>
          dplyr::pull(covar_levs_rast)
        
        # multiply prediction by original raster (correctly matched names)
        tmp_pred_rast_ls[[i]] <- tmp_pred_rast_ls[[i]] * non_bin_r[[lyr_orig_nm]]
        
        names(tmp_pred_rast_ls)[i] <- lyr_nm
        
      }
      
      # %%%%%
      # Update - now that lmer predictions include sd, lo, hi, do below for lmer and brms
      pred_nms <- names(tmp_pred_rast_ls[[1]])
      
      # Sum weighted predictions
      tmp_pred_rast_ls <- purrr::map(.x = pred_nms,
                                     function(.x, 
                                              terra_ls = tmp_pred_rast_ls){
                                       lyr_of_interest <- .x
                                       tmp_terra_ls <- purrr::map(.x = terra_ls,
                                                                  function(.x, 
                                                                           lyr_nm = lyr_of_interest){
                                                                    return(.x[[lyr_nm]])
                                                                  })
                                       tmp_terra <- tmp_terra_ls |>
                                         unlist() |>
                                         c() |>
                                         terra::rast() |>
                                         sum()
                                       
                                       names(tmp_terra) <- lyr_of_interest                    
                                       
                                       return(tmp_terra)
                                     })
      
      data_pred <- tmp_pred_rast_ls |>
        unlist() |>
        c() |>
        terra::rast()
      # %%%%%
      
    }
    
    # otherwise, just combine all into one raster stack
    if (n_prop_covar_lyrs == 0){
      
      tmp_pred_rast <- c(tmp_pred_rast, tmp_bin_rast) #[[1]]
      if (is.list(tmp_pred_rast)){
          tmp_pred_rast <- tmp_pred_rast[[1]]
        }

        assertthat::assert_that(
          inherits(tmp_pred_rast, "SpatRaster")
        )
      
      # Check all covars in new SpatRaster
      assertthat::assert_that(all(fix_vars %in% names(tmp_pred_rast)))
      
      check_covar_levels(obj$model, tmp_pred_rast, fix_vars)
      
      # ranef covar check
      if ( (inherits(proj_specs$re_form, "formula") | 
            is.null(proj_specs$re_form)) ){
        check_covar_levels(obj$model, tmp_pred_rast, re_form_vars)
      }
      
      if (inherits(obj, "HibModelLME4")){
        data_pred <- terra_proj_lme4_pred(obj$model, tmp_pred_rast, backtransform,
                                          proj_specs)
      }
      if (inherits(obj, "HibModelBRMS")){
        data_pred <- terra_proj_brms_pred(obj$model, tmp_pred_rast, backtransform,
                                          proj_specs)
      }
    }
  }
  
  # Optional rescaling
  data_pred_rescaled <- NULL
  if (rescale_pred){
    # divide predictions by baseline val
    data_pred_rescaled <- data_pred / base_val
    
    # Drop some brms estimates of error - doesn't make sens to scale these by baseline...
    # Or does it??
    names(data_pred_rescaled) <- paste0("rescaled_", names(data_pred_rescaled))
  }
  

  # if time in data (original), ensure it is in output
  # Assume that the time for a data SpatRaster is the same for all layers - i.e a year
  t_orig <- terra::time(data)[1]

  data <- c(data, data_pred, data_pred_rescaled)
    
  if (!is.na(t_orig)){
    terra::time(data) <- rep(t_orig, terra::nlyr(data))  
  }
  
  # return SpatRaster
  return(data)
}


#' terra_covar_prep: A helper function to prepare SpatRaster data for projection
#' 
#' @description 
#' This function prepares the supplied SpatRaster to be ready for model projection
#' @param covar_nm A `character`, the name of the model covariate 
#' @param covar_df A `data.frame`, with info on covariate, levels, and the associated
#' names present in the raster data
#' @param data_rast A[`SpatRaster` that can be used to project a model.
#' @param mod_data A `data.frame`, or list (`INLA`), of the data used to fit the 
#' model
#' @param input_type A `character`, either "binary" or "non-binary", indicating 
#' the type of data in the raster - i.e., is each 'grid cell' only a single level of 
#' a covariate, or can it contain proportions of multiple levels.
#' 
#' @return SpatRaster
#' 
#' @keywords data, preparation, model, projection
#' 
#' @seealso 
#' \code{\link{hib_project_model}}
#' 
#' @examples 
#' \dontrun{
#'   r <- terra_covar_prep(covar_nm, covar_df, data_rast, mod_data, input_type = "binary")
#' } 
#' @noRd
terra_covar_prep <- function(covar_nm, 
                             covar_df, 
                             data_rast, 
                             mod_data, 
                             input_type = "binary"){
  # TODO: better term than binary
  
  # Assertions
  assertthat::assert_that(
    is.character(covar_nm) || is.factor(covar_nm), 
    is.data.frame(covar_df), 
    inherits(data_rast, "SpatRaster"), 
    is.data.frame(mod_data) || is.list(mod_data),
    is.character(input_type)
  )
  
  input_type <- match.arg(input_type, c("binary", "non-binary"), several.ok = FALSE)
  
  # Subset to covar of interest
  df <- covar_df |>
    dplyr::filter(covar_nm == covar_nm &
                    !is.na(covar_levs_rast))
  
  r <- data_rast[[df$covar_levs_rast]]
  
  if (input_type == "binary"){
    if (max(terra::values(sum(r)), na.rm = T) != 1){
      message(crayon::yellow("\nWarning: Binary SpatRaster contains cells with >1 category! Layers will be discarded ..."))
      return (NULL)
    }  
  }
  
  # get variable level IDs
  # What if mod_data[[covar_nm]] is not a factor!!?? - manually set... 
  # - or, should this be dealt with previously! - i.e in hib_fit... - for brms
  fac_lev_df <- data.frame(ID = as.numeric(mod_data[[covar_nm]]),
                           category = mod_data[[covar_nm]]) |>
    unique()
  
  
  r <- lapply(r, function(lyr, covar_df = df, lev_df = fac_lev_df){
    
    # Need to convert values > 0 to the number that  matches that factor level 
    # in the model df; e.g. primin = 1 in model levels, needs to be 1 here too 
    covar_level <- covar_df |>
      dplyr::filter(covar_levs_rast == names(lyr)) |>
      dplyr::pull(covar_levs)
    
    covar_level_numeric <- lev_df$ID[lev_df$category == covar_level]
    
    if (input_type == "binary"){
      terra::values(lyr)[terra::values(lyr) > 0] <- covar_level_numeric
    }
    if (input_type == "non-binary"){
      terra::values(lyr) <- covar_level_numeric
    }
    
    # terra::values(lyr)[terra::values(lyr) > 0] <- covar_level_numeric
    # terra::values(lyr) <- covar_level_numeric
    # Or, set all cells to that lyr value - and then rely on subsequent multiplication
    # to remove 0s and NAs...

    if (input_type == "non-binary"){
      lyr <- lyr |>
        terra::as.factor() |>
        terra::categories(value = lev_df)
      # %%%
    }
    return(lyr)
  })
  
  r <- unlist(r) |>
    terra::rast() 
  
  if (input_type == "binary"){

    r <- r |>
      max() |>
      terra::as.factor() |>
      terra::categories(value = fac_lev_df)
    
  }
  
  # %%%
  names(r) <- rep(covar_nm, dim(r)[3]) 
  
  return(r)
}


#' terra_proj_lme4_pred: A helper function to project a model using SpatRaster data
#' 
#' @description 
#' This function projects/predicts the supplied model using the SpatRaster
#'  
#' @param mod An `lme4::merMod` object 
#' @param data A `SpatRaster`, prepared for projection
#' @param backtransform A `function` (or `NULL`) that can backtransform the 
#' response variable to its original/natual scale.
#' @param lme4_specs A named `list` (or `NULL`) containing details relevant to 
#' projecting a `HibModelLME4`. \cr
#' Named items should consist of:
#' - `bootstrap` (`logical`), indicating whether bootstrapped prediction intervals 
#' should be estimated. Can be slow for large models.
#' - `n_boot` (`numeric`), how many bootstrap samples should be used to estimate
#' prediction uncertainty.
#' - `calc_se` (`logical`), specify whether or not a standard error based uncertainty 
#' estimate (via [`marginaleffects::predictions`]) is calculated and used to generate prediction 
#' intervals.
#' - `conf_level` (`numeric`) is used to specify the confidence levels of intervals 
#' generated via bootstrapping or the SE-based approach. Values should fall between
#' 0 and 1, typically in the range 0.5 to 0.95.
#' Multiple values can be provided.
#' - `re_form` (`formula`, `NULL`, or `NA`), as `re.form` in [`predict.merMod()`], 
#' used to specify random effects to condition on when predicting.
#' 
#' @return SpatRaster
#' 
#' @keywords model, projection
#' 
#' @seealso 
#' \code{\link{hib_project_model}}
#' 
#' @examples 
#' \dontrun{
#'   pred <- terra_proj_lme4_pred(mod, data, backtransform, lme4_specs)
#' } 
#' @noRd
terra_proj_lme4_pred <- function(mod, 
                                 data, 
                                 backtransform, 
                                 lme4_specs){
              
  # TODO: silent/verbose arg - use if multiple layers in input raster, or a list...

  # Assertions
  assertthat::assert_that(
    inherits(mod, "merMod"),
    inherits(data, "SpatRaster"), 
    is.null(backtransform) || is.function(backtransform),
    is.list(lme4_specs) || is.null(lme4_specs)
  )
  
  # %%%
  # Project
  # TODO: could also incl df.. - more general
  
  # %%%
  message(crayon::green("\nProjecting model..."))
  data_pred <- terra::predict(data,
                              mod,
                              type = "response", # lme4_specs$type, ??
                              re.form = lme4_specs$re_form, 
                              na.rm = T)
  names(data_pred) <- "prediction_mean"
  
  
  # %%%%%
  if (lme4_specs$bootstrap == TRUE){
    # %%%%
    message(crayon::green(paste("\nGenerating bootstrapped prediciton intervals using",
                                lme4_specs$n_boot, "bootstrap samples...")))
    
    
    tmp_rast_df <- terra::as.data.frame(data, xy = T)

    boot_obj <- lme4::bootMer(mod, 
                          function(m, rast = data){
                            tmp <- terra::predict(rast, m, 
                                                  type = "response", # lme4_specs$type, ??
                                                  re.form = lme4_specs$re_form, 
                                                  na.rm = T) |>
                              terra::as.data.frame()
                            tmp[,1]
                          }, 
                          nsim = lme4_specs$n_boot,
                          .progress = "txt")
    
    # %%%
    # Use conf_level/conf_quants to get prediction_q...
    conf_quants <- conf_level_to_quantiles(lme4_specs$conf_level)
    
    # add 0.5 to conf_quants, unique values
    conf_quants_ <- unique(c(0.5, conf_quants))
    
    # apply, rename - prediction_boot_q..., prediction_boot_median
    uncert_df <- t(apply(boot_obj$t, 2, function(x) as.numeric(stats::quantile(x, probs = conf_quants_, na.rm = TRUE)))) |> 
      data.frame()
    
    colnames(uncert_df) <- paste0("prediction_boot_q", conf_quants_*100)
    uncert_df <- uncert_df |>
      dplyr::rename(prediction_boot_median = prediction_boot_q50)
    
    uncert_df$prediction_boot_sd <- apply(boot_obj$t, 2, stats::sd)
    # %%%

    uncert_df <- dplyr::bind_cols(tmp_rast_df[,c("x", "y")], uncert_df)
    
    data_pred_uncert <- terra::rast(uncert_df,
                                    crs = terra::crs(data))
    
    data_pred <- c(data_pred, data_pred_uncert)

    rm(tmp_rast_df)
    # %%%%
  }
  
  # se-based CIs... and conf level...
  if (lme4_specs$calc_se == T){
    message(crayon::green("\nCalculating standard error of model predictions..."))
    # get conf level if specified - only one, or multiple
    
    # Need to convert data to df...
    tmp_rast_df <- terra::as.data.frame(data, xy = T)
    # if only 1 conf_level - get single df for se CI
    # prediction_lo_<conf_level>
    se_df <- marginaleffects::predictions(mod,
                                          tmp_rast_df,
                                          re.form = lme4_specs$re_form,
                                          type = "response",
                                          conf_level = lme4_specs$conf_level[1]) |> 
      data.frame() |> 
      dplyr::select(estimate, std.error, conf.low, conf.high)
    
    colnames(se_df) <- c("estimate",
                         "prediction_se",
                         paste0("prediction_se_", lme4_specs$conf_level[1], "_lo"),
                         paste0("prediction_se_", lme4_specs$conf_level[1], "_hi"))
    
    # if multiple conf_level - get df once, then use conf level to adjust CI mult and generate se CIs
    if (length(lme4_specs$conf_level) > 1){
      for (i in 2:length(lme4_specs$conf_level)){
        tmp_alpha <- 1 - lme4_specs$conf_level[i]
        tmp_mult <- abs(stats::qnorm(tmp_alpha / 2))
        
        tmp_lo_nm <- paste0("prediction_se_", lme4_specs$conf_level[i], "_lo")
        tmp_hi_nm <- paste0("prediction_se_", lme4_specs$conf_level[i], "_hi")
        
        se_df[, tmp_lo_nm] <- se_df$estimate - (tmp_mult * se_df$prediction_se)
        se_df[, tmp_hi_nm] <- se_df$estimate + (tmp_mult * se_df$prediction_se)
      }
    }
    
    se_df <- se_df |>
      dplyr::select(-estimate)
    
    se_df <- dplyr::bind_cols(tmp_rast_df[,c("x", "y")], se_df)
    
    data_pred_se <- terra::rast(se_df,crs = terra::crs(data))
    
    data_pred <- c(data_pred, data_pred_se)

    rm(tmp_rast_df, se_df)
  }
  
  # %%%
  if (!is.null(backtransform)){
    p_gt_0_orig <- all(terra::values(data_pred, na.rm = T) > 0)
    
    data_pred <- backtransform(data_pred)
    p_gt_0_btr <- all(terra::values(data_pred, na.rm = T) > 0)
    
    if ((p_gt_0_orig & !p_gt_0_btr)){
      message(crayon::yellow("\nWarning: backtransformation of predictions introduced negative values!"))
    }
  }
  
  return(data_pred)
}



# %%%%
#' terra_proj_brms_pred: A helper function to project a model using SpatRaster data
#' 
#' @description 
#' This function projects/predicts the supplied model using the SpatRaster 
#'
#' @param mod A `brms::brmsfit` object 
#' @param data A `SpatRaster`, prepared for projection
#' @param backtransform A `function` (or `NULL`) that can backtransform the 
#' response variable to its original/natual scale.
#' @param brms_specs A named `list` (or `NULL`) containing details relevant to 
#' projecting a `HibModelBRMS`.
#' Named items should consist of:
#' - `incl_obs_uncert` (`logical`), to determine whether or not observation-level 
#' uncertainty is included in predictions.
#' Setting this parameter to `TRUE` will produce predictions using
#' [`predict.brmsfit()`] (larger uncertainty), setting to `FALSE` calls 
#' [`fitted.brmsfit()`] internally.
#' - `conf_level` (`numeric`) is used to specify the confidence levels of intervals 
#' generated. Values should fall between 0 and 1, typically in the range 0.5 to 0.95. 
#' Internally converted to quantiles and passed to the `probs` argument of `predict.brmsfit()`/
#' `fitted.brmsfit()`.
#' - `re_form` (`formula`, `NULL`, or `NA`), as `re_formula` in `predict.brmsfit()`, 
#' used to specify random effects to condition on when predicting.
#' 
#' @return SpatRaster
#' 
#' @keywords model, projection
#' 
#' @seealso 
#' \code{\link{hib_project_model}}
#' 
#' @examples 
#' \dontrun{
#'   pred <- terra_proj_brms_pred(mod, data, backtransform, brms_specs)
#' } 
#' @noRd
terra_proj_brms_pred <- function(mod, 
                                 data, 
                                 backtransform, 
                                 brms_specs){
  
  # TODO: silent/verbose arg - use if multiple layers in input raster, or a list...

  # Assertions
  assertthat::assert_that(
    inherits(mod, "brmsfit"),
    inherits(data, "SpatRaster"), 
    is.null(backtransform) || is.function(backtransform),
    is.list(brms_specs) || is.null(brms_specs)
  )
 
  # %%%
  # Project
  # TODO: could also incl df.. - more general

  # %%%%%
  # %%%%%

  # Convert conf_level to conf_quants
  conf_quants <- conf_level_to_quantiles(brms_specs$conf_level)

  message(crayon::green("\nProjecting model..."))
  if (brms_specs$incl_obs_uncert == T){
    mean_pred <- terra::predict(data,
                                mod,
                                re_form = brms_specs$re_form,
                                probs = conf_quants,
                                na.rm = T)
    med_pred <- terra::predict(data,
                               mod,
                               re_form = brms_specs$re_form,
                               robust = T,
                               na.rm = T)
    
    names(mean_pred)[names(mean_pred) == "Estimate"] <- "prediction_mean"
    names(mean_pred)[names(mean_pred) == "Est.Error"] <- "prediction_sd"

    names(mean_pred) <- gsub("Q", "prediction_q", names(mean_pred), fixed = T)
  
    names(med_pred)[names(med_pred) == "Estimate"] <- "prediction_median"
    names(med_pred)[names(med_pred) == "Est.Error"] <- "prediction_mad"
    
    data_pred <- c(mean_pred, med_pred[[c("prediction_median", "prediction_mad")]])
    
    if (!is.null(backtransform)){
      data_pred <- backtransform(data_pred)
    }
  }
  
  # %%%%%
  # And for incl_obs_uncert == F
  # Needs fitted, uses a df...
  if (brms_specs$incl_obs_uncert == F){
    tmp_df <- terra::as.data.frame(data, xy = T)
    mean_pred <- stats::fitted(mod, tmp_df, 
                        scale = "response", # brms_specs$scale
                        re_form = brms_specs$re_form,
                        probs = conf_quants) |>
      as.data.frame()

    median_pred <- stats::fitted(mod, tmp_df, 
                          scale = "response", # brms_specs$scale
                          re_form = brms_specs$re_form,
                          robust = T) |>
      as.data.frame()
    
    
    mean_pred <- mean_pred |>
      dplyr::rename(prediction_mean  = Estimate,
                    prediction_sd = Est.Error) |>
      dplyr::rename_with(~gsub("Q", "prediction_q", .x, fixed = T))

    median_pred <- median_pred |>
      dplyr::rename(prediction_median  = Estimate,
                    prediction_mad = Est.Error)
    tmp_pred <- dplyr::bind_cols(tmp_df[,c("x", "y")], 
                                 mean_pred,
                                 median_pred |>
                                   dplyr::select(prediction_median, 
                                                 prediction_mad))
    
    data_pred <- terra::rast(tmp_pred,
                             crs = terra::crs(data))
    
    if (!is.null(backtransform)){
      p_gt_0_orig <- all(terra::values(data_pred, na.rm = T) > 0)

      data_pred <- backtransform(data_pred)

      p_gt_0_btr <- all(terra::values(data_pred, na.rm = T) > 0)

      if ((p_gt_0_orig & !p_gt_0_btr)){
        message(crayon::yellow("\nWarning: backtransformation of predictions introduced negative values!"))
      }
    }
  }
  
  return(data_pred)
}

# TODO:
# terra_proj_inla_pred <- function(obj, data, backtransform, rescale_pred){
#   # Assertions...
#   
# }   


#' terra_to_df_inla: A helper function to wrap the preparation of SpatRaster data 
#' and then return a dataframe - for use in INLA model projections
#' 
#' @description 
#' This function prepares a SpatRaster for use in INLA projections
#' 
#' @param data A `SpatRaster` that can be used to project a model.
#' @param obj A `HibModelINLA` object.
#' @param fix_vars A `character` vector of model fixed effects.
#' @param rand_vars A `character` vector of model random effects.
#' @param re_form_vars A `character` of model random effects to condition upon 
#' when predicting.
#' 
#' @return data.frame
#' 
#' @keywords model, projections
#' 
#' @examples 
#' \dontrun{
#'   tmp_df <- terra_to_df_inla(obj, data, fix_vars, rand_vars, re_form_vars)
#' } 
#' @noRd
terra_to_df_inla <- function(data, obj, fix_vars, rand_vars, re_form_vars){

  # TODO: maybe co-opt to use fitted for brms rather then predict - for terra input
  
  # assertions
  assertthat::assert_that(
    inherits(obj, "HibModelINLA"),
    inherits(data, "SpatRaster"), 
    is.character(fix_vars),
    is.character(rand_vars) || is.na(c(rand_vars)) || is.null(rand_vars),
    is.character(re_form_vars) || is.na(c(re_form_vars)) || is.null(re_form_vars)
  )
  # Covar checks
  
  if (all(fix_vars %in% names(data))){
    # also check for all ranefs... re_form_vars
    if (all(re_form_vars %in% names(data))){
      # convert terra to df...
      data_df <- terra::as.data.frame(data, xy = T)
      } else {
        if (is.na(re_form_vars)){
          data_df <- terra::as.data.frame(data, xy = T)
        }
        if (!is.na(re_form_vars)){
          # Get missing re_form_vars
          missing_covars <- re_form_vars[!re_form_vars %in% names(data)]
          missing_covars_chr <- paste(missing_covars, collapse = " ; ")

          message(crayon::red(paste0("\nError: random effect variable specified in ",
                                    "re_form is missing.\n")))
          message(crayon::red(paste("Missing covariates: ",
                                missing_covars_chr, sep = "")))
          stop()  
        }
        
      }
    
  } else {
    # Any layer with a name in fix_vars is good to use 
    covars_in_rast <- names(data)[names(data) %in% fix_vars]
    # also get ranefs. - assume not proportions/binary - as in model
    if (!all(is.na(re_form_vars))){
      ranefs_in_rast <- names(data)[names(data) %in% re_form_vars]
      covars_in_rast <- c(covars_in_rast, ranefs_in_rast)
    }

    tmp_pred_rast <- NULL
    n_pred_rast_lyrs <- 0
    if (length(covars_in_rast) > 0 ){
      tmp_pred_rast <- data[[covars_in_rast]]  
      n_pred_rast_lyrs <- dim(tmp_pred_rast)[3]
    }
    
    # any layer with a name not in fix_covars needs to be aligned to a var
    non_covars_in_rast <- names(data)[!names(data) %in% fix_vars]
    # drop ranefs in data from non_covars_in_rast
    # Easy way to drop ranefs_in_rast 
    non_covars_in_rast <- non_covars_in_rast[!non_covars_in_rast %in% covars_in_rast]
    # But, also want to drop ranefs that aren't in rast from non_covars_in_rast
    non_covars_in_rast <- non_covars_in_rast[!non_covars_in_rast %in% rand_vars]

    tmp_rast <- NULL
    if (length(non_covars_in_rast) > 0 ){
      tmp_rast <- data[[non_covars_in_rast]]
    }
    
    # Assume that layer names match levels of covariates
    # e.g "LUIPrimary_minimal use" or just "Primary_minimal use"
    covars_not_in_rast <- fix_vars[!fix_vars %in% names(data)]
    
    if (length(non_covars_in_rast) > 0 ){
      # for each covar, get unique levels (if factor in mod_df)

      # Use model-fitting data from obj$model...
      # inla: $.args$data
      
      # inla fit data,
      if (inherits(obj$model, "inla")){
        fit_dat <- obj$model$.args$data
      }
      
      covars_not_in_rast <- purrr::map_dfr(covars_not_in_rast,
                                           function(x){
                                             if (is.factor(fit_dat[[x]])){
                                               tmp_out <- data.frame(covar_nm = x,
                                                                     covar_levs = levels(fit_dat[[x]]))
                                               # Could add covar_lev_idx here... - save doing later...
                                             }
                                             else {
                                               tmp_out <- data.frame(covar_nm = x,
                                                                     covar_levs = NA)
                                             }
                                             tmp_out
                                           })
      
      covars_not_in_rast <- covars_not_in_rast |>
        # Combine covar nm with levels
        dplyr::mutate(covar_nm_levs = paste0(covar_nm, covar_levs), 
                      # check for presence in tmp_rast
                      covar_levs_in_rast = covar_levs %in% names(tmp_rast), 
                      covar_nm_levs_in_rast = covar_nm_levs %in% names(tmp_rast),
                      # Add level name in raster
                      covar_levs_rast = dplyr::case_when(covar_levs_in_rast == T ~ covar_levs,
                                                         covar_nm_levs_in_rast == T ~ covar_nm_levs,
                                                         T ~ NA))
      
      # Check for binary layers...
      is_binary_layer <- purrr::map_dfr(names(tmp_rast),
                                        function(nm){
                                          data.frame(covar_levs_rast = nm,
                                                     is_binary_layer = all(unique(terra::values(tmp_rast[[nm]])) %in% c(0,1, NA, NaN))
                                          )
                                        })
      
      covars_not_in_rast <- covars_not_in_rast |>
        dplyr::left_join(is_binary_layer, by = "covar_levs_rast") |>
        dplyr::group_by(covar_nm) |>
        dplyr::mutate(is_binary_covar = all(is_binary_layer, na.rm = T))
      
      # Split binary and non binary layers
      bin_covars <- covars_not_in_rast |>
        dplyr::filter(is_binary_covar == T)
      
      non_bin_covars <- covars_not_in_rast |>
        dplyr::filter(is_binary_covar == F)
      
      if (length(unique(non_bin_covars$covar_nm)) > 1){
        message(crayon::yellow("\nWarning: SpatRaster with multiple non-binary layers detected! Projections are not yet supported...\n"))
        stop("In development")
      }
      
      # For binary layers - factor is layer name, minus covar nm if name is covar_nm_lev
      # Collapse to single stack with layer names as covar_nm
      
      # only if binary layers
      tmp_bin_rast <- NULL
      n_bin_covar_layers <- 0
      if (nrow(bin_covars) > 0){
        tmp_bin_rast <- purrr::map(unique(bin_covars$covar_nm),
                                   terra_covar_prep,
                                   covar_df = bin_covars,
                                   data_rast = tmp_rast,
                                   mod_data = fit_dat,
                                   input_type = "binary") |>
          unlist() |>
          c() |>
          terra::rast()
        
        n_bin_covar_layers <- dim(tmp_bin_rast)[3]
      }
      
      # For non-binary layers - factor is layer name, minus covar nm if name is covar_nm_lev
      # BUT, also need to record values for later aggregation
      # currently assume only one non_bin_covar...
      
      tmp_prop_rast <- NULL
      n_prop_covar_lyrs <- 0
      if (nrow(non_bin_covars) > 0){
        # No need for purrr - only gets here if 1 non-binary covar...
        tmp_prop_rast <- purrr::map(unique(non_bin_covars$covar_nm),
                                    terra_covar_prep,
                                    covar_df = non_bin_covars,
                                    data_rast = tmp_rast,
                                    mod_data = fit_dat,
                                    input_type = "non-binary") |>
          unlist() |>
          c() |>
          terra::rast()
        # Want to end up with N raster stacks - 1 per prop level
        n_prop_covar_lyrs <- dim(tmp_prop_rast)[3]
      }
      # Check tmp_prop_rast is a terra SpatRaster
      
    }
    
    # check N layers of each
    # if n_prop_covar_lyrs > 1  - list of rasters - each named based on covar level
    if (n_prop_covar_lyrs > 1){
      
      non_bin_sub <- non_bin_covars |>
        dplyr::filter(!is.na(covar_levs_rast))
      
      non_bin_r <- tmp_rast[[non_bin_sub$covar_levs_rast]]
      
      tmp_pred_df_ls <- vector("list", n_prop_covar_lyrs)
      names(tmp_pred_df_ls) <- rep("", n_prop_covar_lyrs)
      for (i in 1:n_prop_covar_lyrs){
        # need to be careful combining these, if first item is null,
        # output becomes a list of 1...

        # %%%
        tmp_ <- c(tmp_pred_rast, tmp_bin_rast, tmp_prop_rast[[i]])
        # If tmp is a list, get [[1]]
        if (is.list(tmp_)){
          tmp_ <- tmp_[[1]]
        }

        assertthat::assert_that(
          inherits(tmp_, "SpatRaster")
        )
        # If is a SpatRaser, use as is.
        tmp_r <- tmp_ 
        
        # %%%%
        # Check fix vars and re-form_vars are present
        assertthat::assert_that(all(fix_vars %in% names(tmp_r)))
        check_covar_levels(obj$model, tmp_r, fix_vars)
        
        if (!is.na(re_form_vars)){
          check_covar_levels(obj$model, tmp_r, re_form_vars)
        }        
        # %%%%

        # %%%
        # Convert to df
        tmp_df <- terra::as.data.frame(tmp_r,
                                       xy = T)
        
        # Get proportional weighting
        lyr_nm <- terra::cats(tmp_prop_rast)[[i]] |> 
          dplyr::filter(ID == as.vector(unname(unique(terra::values(tmp_prop_rast[[i]])))) ) |> 
          dplyr::pull(category)
        
        lyr_orig_nm <- non_bin_sub |>
          dplyr::filter(covar_levs == lyr_nm) |>
          dplyr::pull(covar_levs_rast)
        
        prop_df <- terra::as.data.frame(non_bin_r[[lyr_orig_nm]],
                                        xy = T)
        
        prop_df$prop_cover <- prop_df[[lyr_orig_nm]]
        # rename lyr_orig_nm as prop_cover
        tmp_pred_df_ls[[i]] <- dplyr::left_join(tmp_df, 
                                                prop_df[,c("x", "y", "prop_cover")],
                                                by = c("x", "y"))
        
      }
      
      data_df <- dplyr::bind_rows(tmp_pred_df_ls)
      
    }
    
    # otherwise, just combine all into one raster stack
    if (n_prop_covar_lyrs == 0){
      
      tmp_r <- c(tmp_pred_rast, tmp_bin_rast)[[1]]
      # Convert to df
      data_df <- terra::as.data.frame(tmp_r,
                                      xy = T)
    }
  }
  
  return(data_df)  
  
}


#' check_covar_levels: An internal helper function to check that the data provided for the projection 
#' does not contain any factor-levels that are not in the model fitting data.
#' 
#' @description 
#' This function checks whether the supplied projection data is compatible with
#' the supplied model.
#' Returns either nothing (NULL, input data passes checks), or throws an error 
#' (input data fails checks)
#' 
#' @param mod A fitted model, either from `lme4`, `INLA` or `brms`.
#' @param proj_dat The prepared data used to predict from
#' @param covars The fixed-effect covariates in the model # TODO: could also be ranefs
#' 
#' @return NULL
#' 
#' @keywords model projections, data checks
#' 
#' @seealso 
#' \code{\link{hib_project_model}}
#' 
#' @examples 
#' \dontrun{
#'   check_covar_levels(model, projection_df, fixed_covars)
#' } 
#' @noRd
check_covar_levels <- function(mod, proj_dat, covars){
  
  # Assertion checks
  assertthat::assert_that(
    inherits(mod, "lmerMod") || inherits(mod, "glmerMod") ||
      inherits(mod, "brmsfit") || inherits(mod, "inla"),
    is.data.frame(proj_dat) || inherits(proj_dat, "SpatRaster"),
    is.character(covars)
  )
  
  # First check - that all covars are in proj_dat
  if (!all(covars %in% names(proj_dat))){
    
    # Get missing covars
    missing_covars <- covars[!covars %in% names(proj_dat)]
    missing_covars_chr <- paste(missing_covars, collapse = " ; ")
    message(crayon::red("\nError: covariates missing from projection data!"))
    message(crayon::red(paste("\nMissing covariates: ",
                              missing_covars_chr, sep = "")))
    
    stop()
  }

  # store info on data type
  proj_dat_type <- NA
  if (is.data.frame(proj_dat)){
    proj_dat_type <- "df"
  }
  if (inherits(proj_dat, "SpatRaster")){
    proj_dat_type <- "rast"
  }
  
  # lme4 fit data,
  if ((inherits(mod, "lmerMod") || inherits(mod, "glmerMod"))){
    fit_dat <- mod@frame
  }
  # brms fit data
  if (inherits(mod, "brmsfit")){
    fit_dat <- mod$data
  }
  # inla fit data
  if (inherits(mod, "inla")){
    fit_dat <- mod$.args$data
  }
  
  problem_vars <- list()
  
  is_problem_var <- FALSE
  for (i in 1:length(covars)){
    # Check if var is a factor/chr
    var_i <- covars[i]
    if ( (is.factor(fit_dat[[var_i]]) | is.character(fit_dat[[var_i]])) ){
      # Get unique chr values for fit data
      var_i_fit_levs <- fit_dat[[var_i]] |> 
        unique() |>
        as.character()
      # Alternative - get levels directly - if factor
      if (is.factor(fit_dat[[var_i]])){
        var_i_fit_levs <- levels(fit_dat[[var_i]])
      }
      
      # Get unique chr values for proj data 
      var_i_proj <- proj_dat[[var_i]]
      
      if (proj_dat_type == "df"){
        var_i_proj_levs <- var_i_proj |> 
          unique() |>
          as.character()
        # Alternative - get levels directly - if factor
        if (is.factor(var_i_proj)){
          var_i_proj_levs <- levels(var_i_proj)
        }
      }
      if (proj_dat_type == "rast"){
        var_i_proj_levs <- terra::levels(var_i_proj)[[1]][["category"]]
      }

      # Check all var_i_proj_levs are in var_i_fit_levs
      if (!all(var_i_proj_levs %in% var_i_fit_levs)){
        # if not, add to problem_vars
        is_problem_var <- TRUE
        
        problem_vars[[i]] <- dplyr::tibble(Variable = var_i,
                                           Problem_levels = list(var_i_proj_levs[!var_i_proj_levs %in% var_i_fit_levs]) )
      }
      
    }
  }
  
  if (!is_problem_var){
    return(NULL)
  }
  if (is_problem_var){
    # Print an error message - detail variables and levels that are issues...
    message(crayon::red(paste0("\nSome categorical covariates in the projection data ",
                               "are not in the model fitting data!\n\n",
                               "Please fix this, see details below:\n"
                               ) ))
    problem_vars |>
      dplyr::bind_rows() |>
      insight::format_table(use_symbols = T, digits = 0) |>
      insight::export_table() |>
      message()
    
    stop("Unknown covariate levels, see above.", call. = F)
  }
}

# compare data to use for projection with model fitting data - unique values of 
# categorical covars
# input1: fited model object (just model, not hib model object) (or just model fitting df)
# input2: projection df or projection raster/terra
# fix_vars: need to no covariates of interest - 
# access model@frame
# check for factor columns - 
# for each, get levels
# then check the levels of the projection data - 

# If non-compatible levels - drop? from df
# Can they also be dropped from terra rast?
# Maybe better to throw an error - not possible to make predictions...
# print non-compatible levels

# Also need to check for NA values - should drop these I think... 
# or throw a warning/erorr



#' conf_level_to_quantiles: An internal helper function to convert confidence 
#' levels to associated quantiles.
#' 
#' @description 
#' This function converts confidence levels to the associated quantiles. e.g., 
#' converting c(0.5, 0.9) to c(0.05 0.25 0.75 0.95)
#' 
#' @details 
#' Expects a numeric vector
#' 
#' @param x A `numeric` vector of confidence levels between 0 and 1, typically
#' 0.5 to 0.95.
#' 
#' @return numeric
#' 
#' @keywords model projections, data checks
#' 
#' @seealso 
#' \code{\link{hib_project_model}}
#' 
#' @examples 
#' \dontrun{
#'   conf_level_to_quantiles(c(0.5, 0.9, 0.95))
#' } 
#' @noRd
conf_level_to_quantiles <- function(x){
  assertthat::assert_that(
    is.numeric(x),
    min(x) > 0,
    max(x) < 1
  )
  
  lwr <- 0.5 - (x/2)
  upr <- 0.5 + (x/2)
  
  return(c(sort(lwr), upr))
}


#' check_proj_specs: An internal helper function to check model projection specs.
#' 
#' @description 
#' This function checks and corrects model projection specs.
#' If specifications are missing, default options are added.
#' 
#' @param specs A named `list` of model specs - for details see `lme4_specs`, 
#' `inla_specs` and `brms_specs` arguments to `hib_project_model.`
#' @param mod_type A `character` indicating the model type.
#' 
#' @return list
#' 
#' @keywords model projections, data checks
#' 
#' @seealso 
#' \code{\link{hib_project_model}}
#'  
#' @examples 
#' \dontrun{
#'   check_proj_specs(lme4_specs, "lme4")
#' } 
#' @noRd
check_proj_specs <- function(specs, mod_type){
  
  assertthat::assert_that(
    is.list(specs) || is.null(specs),
    is.character(mod_type)
  )
  # Specify default specs
  if (mod_type %in% c("HibModelLME4", "lme4")){
    default_specs <- list(bootstrap  = FALSE, 
                          n_boot     = NA, 
                          re_form    = NA, 
                          calc_se    = TRUE, 
                          conf_level = 0.95)
    mod_type_nm <- "lme4"
  }
  if (mod_type %in% c("HibModelINLA", "inla")){
    default_specs <- list(re_form      = NA, 
                          link_inverse = NULL)
    mod_type_nm <- "inla"
  }
  if (mod_type %in% c("HibModelBRMS", "brms")){
    default_specs <- list(incl_obs_uncert = FALSE,
                          re_form         = NA, 
                          conf_level      = 0.95)
    mod_type_nm <- "brms"
  }
  if (mod_type == "rescale"){
    default_specs <- list(rescale          = FALSE, 
                          baseline_covar   = NULL, 
                          baseline_re_form = NA) 
    mod_type_nm <- "rescale"
  }
  
  # If specs is.null(), set to defaults
  if (is.null(specs)){
    message(crayon::yellow(paste0("\nWarning: ", mod_type_nm, "_specs missing, using default values:\n",
                                  "\t ",
                                  paste(paste(names(default_specs), 
                                              default_specs, 
                                              sep = " = "),
                                        collapse = ",\n\t "),
                                  " ...")))
    return(default_specs)
  }
  
  
  # Assertions to check data type
  if ("bootstrap" %in% names(specs)){
    assertthat::assert_that(is.logical(specs$bootstrap))
  }
  if ("n_boot" %in% names(specs)){
    assertthat::assert_that(is.numeric(specs$n_boot))
  }
  if ("re_form" %in% names(specs)){
    assertthat::assert_that(
      is.na(c(specs$re_form)) || 
        is.null(specs$re_form) || 
        inherits(specs$re_form, "formula")
    )
  }
  if ("calc_se" %in% names(specs)){
    assertthat::assert_that(is.logical(specs$calc_se))
  }
  if ("conf_level" %in% names(specs)){
    assertthat::assert_that(is.numeric(specs$conf_level))
  }
  if ("link_inverse" %in% names(specs)){
    assertthat::assert_that(
      is.null(specs$link_inverse) || 
        is.function(specs$link_inverse)
    )
  }
  if ("incl_obs_incert" %in% names(specs)){
    assertthat::assert_that(is.logical(specs$incl_obs_incert))
  }
  
  
  if ("rescale" %in% names(specs)){
    assertthat::assert_that(
      is.na(c(specs$rescale)) || 
        is.null(specs$rescale) || 
        is.logical(specs$rescale)
    )
    if (is.na(specs$rescale) | is.null(specs$rescale)){
      specs$rescale <- FALSE
    }
  }
  if ("baseline_covar" %in% names(specs)){
    assertthat::assert_that(
      is.null(specs$baseline_covar) || 
        is.data.frame(specs$baseline_covar)
    )
    if (is.data.frame(specs$baseline_covar)){
      assertthat::assert_that(
        nrow(specs$baseline_covar) == 1
      )
    }
  }
  if ("baseline_re_form" %in% names(specs)){
    assertthat::assert_that(
      is.na(c(specs$baseline_re_form)) || 
        is.null(specs$baseline_re_form) || 
        inherits(specs$baseline_re_form, "formula")
    )
  }
  
  
  # Fill in any missing specs
  for (nm in names(default_specs)){
    if (!nm %in% names(specs)){
      # If bootstrap is TRUE, ensure n_boot is numeric
      if (nm == "n_boot"){
        if (specs$bootstrap){
          message(crayon::yellow(paste0("\nWarning: ", nm, " missing from ", mod_type_nm, "_specs, using default value of ",  50, " ...")))
          specs[[nm]] <- 50  
        }
      } else if (nm == "baseline_covar"){
        message(crayon::yellow("\nWarning: ", nm, " missing from ", mod_type_nm, "_specs, using default value of NULL. \n          Model intercept will be used as the baseline."))
        specs[[nm]] <- NULL
      } else if (nm == "link_inverse"){
        message(crayon::yellow("\nWarning: ", nm, " missing from ", mod_type_nm, "_specs, using default value of NULL."))
        specs[[nm]] <- NULL
      } else {
        message(crayon::yellow(paste0("\nWarning: ", nm, " missing from ", mod_type_nm, "_specs, using default value of ",  default_specs[[nm]], " ...")))
        specs[[nm]] <- default_specs[[nm]]
      }
    }
  }
  
  if (mod_type == "rescale"){
    if (specs$rescale == FALSE){
      rescale_specs$baseline_covar <- NULL
      rescale_specs$baseline_re_form <- NA
    }
  }
  
  return(specs)
}
