#' Project a hib, mixed-effects model
#' 
#' @description 
#' This function projects a fitted hib model using the supplied model object and
#' data to be used in projection - i.e., a dataframe/tibble or raster stack of covariates
#' @details 
#' Expects a [`hib_model`] object.
#' @param mod_obj A [`hib_model`] object - either lme4, inla, or brms.
#' @param data Either a dataframe-like object or a raster stack that can be used
#' to project onto
#' @param backtransform A [`function`] that will backtransform the model response
#' to the scale of the original variable, e.g., exp() as the inverse of log()
#' @param rescale_pred A [`logical`] value - should the predictions be rescaled relative
#' to the model intercept (assumes intercept is for baseline environment/LU)
#' @param lme4_specs A list containing details relevant to projecting an lme4-style
#' model. Specifically, named items should consist of 'bootstrap' [`logical`] and 
#' 'n_boot' [`numeric`], indicating whether bootstrapped prediction intervals should be 
#' estimated, and if so, how many bootstrap samples should be used.
#' Additionally, 're_form', as 're.form' in lme4:::predict.merMod(), used to specify random
#' effects to condition on when predicting an lme4 model (formula, NULL, or NA).
#' @param inla_specs Currently unused, but may be useful in future dev.
#' @param brms_specs A list containing details relevant to projecting a brms-style
#' model. Specifically, named item should consist of 'incl_obs_uncert' [`logical`] 
#' to determine whether or not observation-level uncertainty is included in predictions.
#' Setting this parameter to TRUE will produce predictions with larger uncertainty.
#' TODO: also include re_form - as 're_formula' in brms predict/fitted - use as
#' in lme4.
#' @keywords model, projections
#' @examples 
#' \dontrun{
#'   proj <- hib_mod_project(hib_model, proj_data)
#' } 
#' @export
# TODO: replace obj with mod_obj
hib_mod_project <- function(mod_obj, 
                            data,
                            backtransform = NULL,
                            rescale_pred = FALSE,
                            # incl_obs_uncert = F,
                            # lme4_boot = T,
                            # n_boot = 100
                            # TODO: collate model type arguments into: - fix argument documentation
                            lme4_specs = NULL, #list(bootstrap = F, n_boot = NA, re_form = NA),
                            inla_specs = NULL, #list(re_form = ??),
                            brms_specs = NULL, #list(incl_obs_uncert = F, TODO, re_form = as lme4...)
                            # more scope for adjustment, and neater
                            calc_aoa = FALSE,
                            cv_obj = NULL
                            ){
  # Assertions 
  # is data a df/tibble/ raster/stack 
  assertthat::assert_that(
    inherits(mod_obj, "hib_lme4_model") || inherits(mod_obj, "hib_inla_model") || inherits(mod_obj, "hib_brms_model"), # obj is a hib model
    is.data.frame(data) || inherits(data, "SpatRaster") || ( is.list(data) & !is.data.frame(data) & all(unlist(lapply(data, function(x){inherits(x, "SpatRaster")}))) ), 
    is.function(backtransform) || is.null(backtransform),
    is.logical(rescale_pred),
    # is.logical(incl_obs_uncert),
    # is.logical(lme4_boot),
    # is.numeric(n_boot) || is.null(n_boot)
    is.list(lme4_specs) || is.null(lme4_specs),
    is.list(inla_specs) || is.null(inla_specs),
    is.list(brms_specs) || is.null(brms_specs),
    is.logical(calc_aoa),
    inherits(cv_obj, "hib_cv_res") || is.null(cv_obj)
  )
  

  # %%%%%
  # TODO: - incl calculation of DI/AOA
  # Need to also incl cv res
  if (calc_aoa){
    assertthat::assert_that(inherits(cv_obj, "hib_cv_res"))  
    message(crayon::green("Estimating Area-of-Applicability...\n"))
    aoa_res <- hib_calc_aoa(mod_obj, cv_obj, data)
  }
  if (!calc_aoa){
    aoa_res <- NULL
  }
  # %%%%%


  # Check type of model
  # %%%%%
  # LME4
  # %%%%%
  if (inherits(mod_obj, "hib_lme4_model")){
    message(crayon::green("lme4 model detected! Projection in progress...\n"))
    # TODO: wrap in hib_mod_proj_lme4??
    
    # %%%%%
    # Get info re bootstrapping
    if (!is.null(lme4_specs)){
      # Check for bootstrap - logical
      if ("bootstrap" %in% names(lme4_specs)){
        assertthat::assert_that(is.logical(lme4_specs$bootstrap))
      }
      # if absent, set to F
      if (!"bootstrap" %in% names(lme4_specs)){
        lme4_specs$bootstrap <- FALSE
      }
      
      
      # Check for n_boot - NA/integer/numeric
      if (lme4_specs$bootstrap){
        if ("n_boot" %in% names(lme4_specs)){
          assertthat::assert_that(is.numeric(lme4_specs$n_boot))
          # Could fix n_boot to 50 if is.na(n_boot) but bootstrap == TRUE
        }
        # if absent, set to F
        if (!"n_boot" %in% names(lme4_specs)){
          message(crayon::yellow("Warning: n_boot missing, using default value of 50...\n"))
          lme4_specs$n_boot <- 50
        }
      }
      # if absent, if bootstrap, set to 50, if not bootstrap, set to NA
      if (lme4_specs$bootstrap == FALSE){
        lme4_specs$n_boot <- NA
      }

      if ("re_form" %in% names(lme4_specs)){
        assertthat::assert_that(
          is.na(c(lme4_specs$re_form)) || is.null(lme4_specs$re_form) || inherits(lme4_specs$re_form, "formula")
        )
      }
      if (!"re_form" %in% names(lme4_specs)){
        message(crayon::yellow("Warning: re_form missing, using default value of NA...\n"))
        lme4_specs$re_form <- NA
      }

    }
    if (is.null(lme4_specs)){
      message(crayon::yellow(paste0("Warning: lme4_specs missing, using default values:\n",
                                    "\tbootstrap = FALSE,\n",
                                    "\tn_boot = NA,\n",
                                    "\tre_form = NA",
                                    " ...\n")))
      lme4_specs <- list(bootstrap = FALSE, n_boot = NA, re_form = NA) #
    }
    
    # print(lme4_specs)
    # %%%%%
    
    # Get baseline value for optional rescaling
    
    # %%%%
    # get model intercept
    if ("(Intercept)" %in% mod_obj$fixedef_coefficients$Coefficient_name){
      base_val <- mod_obj$fixedef_coefficients$Coefficient_value[mod_obj$fixedef_coefficients$Coefficient_name == "(Intercept)"]
    } else {
      base_val <- mod_obj$fixedef_coefficients$Coefficient_value[1]
    }
    # %%%%
    # get model intercept
    # base_val <- obj$fixedef_coefficients$Coefficient_value[obj$fixedef_coefficients$Coefficient_name == "(Intercept)"]
    # if backtrans - backrans intercept to get baseline val
    if (!is.null(backtransform)){
      base_val <- backtransform(base_val)
    }
    
    message(crayon::green("Preparing input data for projection...\n"))
    # Get some basic info on model covars
    all_covars <- all.vars(formula(delete.response(terms(formula(mod_obj$model)))))
    # But what about random effects... - 
    rand_vars <- names(lme4::ranef(mod_obj$model))
    fix_vars <- all_covars[!all_covars %in% rand_vars]
    
    # TODO: lme4_spec checks
    # lme4_specs - if null, set to default (F, NA)
    # otherwise, assert that bootstrap is logical and n_boot is numeric
    # then just pass checked lme4_specs to functions - terra_proj...
    # terra_proj_lme4_mod...
    
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
      
      # need to check any random vars in lme4_specs$re_form
      if (inherits(lme4_specs$re_form, "formula")){
        # print("Checking random variables")
        re_form_vars <- all.vars(lme4_specs$re_form)
        assertthat::assert_that(all(re_form_vars %in% rand_vars), 
                                msg = "Specified re_form structure contains variables not in the original model.")
        # all.vars(lme4_specs$re_form[[2]])
        check_covar_levels(mod_obj$model, data, re_form_vars)
      }
      if (is.null(lme4_specs$re_form)){
        check_covar_levels(mod_obj$model, data, rand_vars)
      }

      message(crayon::green("Projecting model...\n"))
      pred_df <- as.data.frame(predict(mod_obj$model,
                                       data,
                                       type = "response", # lme4_specs$type,
                                       re.form = lme4_specs$re_form)) 
      colnames(pred_df) <- "prediction_mean" # prediction
      
      # TODO: current approach to backtransforming/rescaling sd doesn't make sense...

      # %%%
      # Use bootMer to get uncertainty
      if (lme4_specs$bootstrap == T){
        message(crayon::green("Bootstrapping model predictions...\n"))
        boot_obj <- lme4::bootMer(mod_obj$model, 
                                  function(mod, df = data){predict(mod, df, 
                                    type = "response", # lme4_specs$type,
                                    re.form = lme4_specs$re_form)}, 
                                  nsim = lme4_specs$n_boot)
        
        pred_df$prediction_median <- apply(boot_obj$t, 2, function(x) as.numeric(quantile(x, probs = 0.5, na.rm = TRUE)))
        pred_df$prediction_sd <- apply(boot_obj$t, 2, sd)
        pred_df$prediction_q2.5 <- apply(boot_obj$t, 2, function(x) as.numeric(quantile(x, probs = 0.025, na.rm = TRUE)))
        pred_df$prediction_q97.5 <- apply(boot_obj$t, 2, function(x) as.numeric(quantile(x, probs = 0.975, na.rm = TRUE)))
        # std.err <- apply(boot_test$t, 2, sd)
        # CI.lower1 <- pred - std.err*1.96
        # CI.upper1 <- pred + std.err*1.96
        # %%%
        # Should back trans and rescale here...
      }
      
      # What if we want to include a random effect?? - lme4_specs$re_form
      
      if (!is.null(backtransform)){
        pred_df <- backtransform(pred_df)
      }
      pred_rescaled_df <- NULL
      if (rescale_pred){
        # divide predictions by baseline val
        pred_rescaled_df <- pred_df / base_val
        colnames(pred_rescaled_df) <- paste("rescaled", colnames(pred_rescaled_df), sep = "_")
      }
      
      data <- dplyr::bind_cols(data, pred_df, pred_rescaled_df)
      
    }
    
    # Single raster object
    if (inherits(data, "SpatRaster")){
      # print("running terra proj")
      data <- terra_proj(data, mod_obj, backtransform, rescale_pred, base_val, fix_vars,
                          rand_vars,
                         lme4_specs = lme4_specs)
    }
    # List of raster objects
    if ((is.list(data) & !is.data.frame(data))){
      # & inherits(data[[1]], "SpatRaster")
      # Apply individual terra proj to each item of list...
      data <- purrr::map(data,
                         terra_proj,
                         obj = mod_obj, 
                         backtransform = backtransform, 
                         rescale_pred = rescale_pred, 
                         base_val = base_val, 
                         fix_vars = fix_vars,
                         rand_vars = rand_vars,
                         lme4_specs = lme4_specs)
                         # lme4_boot = lme4_boot,
                         # n_boot = n_boot)
    }
    
  } 
  # %%%%%
  # INLA
  # %%%%%
  else if (inherits(mod_obj, "hib_inla_model")){

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
    message(crayon::green("INLA model detected! Projection in progress...\n"))
    # stop("In development")
    
    # %%%
    # Get info re bootstrapping
    if (!is.null(inla_specs)){
    
      if ("re_form" %in% names(inla_specs)){
        assertthat::assert_that(
          is.na(c(inla_specs$re_form)) || is.null(inla_specs$re_form) || inherits(inla_specs$re_form, "formula")
        )
      }
      if (!"re_form" %in% names(inla_specs)){
        message(crayon::yellow("Warning: re_form missing, using default value of NA...\n"))
        inla_specs$re_form <- NA
      }

    }
    if (is.null(inla_specs)){
      message(crayon::yellow(paste0("Warning: inla_specs missing, using default values:\n",
                                    "\tre_form = NA",
                                    " ...\n")))
      inla_specs <- list(re_form = NA) #
    }
    # print(inla_specs)
    # %%%



    # Get baseline value for optional rescaling
    # get model intercept
    int_id <- which(grepl("intercept", rownames(mod_obj$model$summary.fixed), ignore.case = T))
    # TODO: What if another covar also contains intercept...
    if (length(int_id) > 1){
      message(crayon::yellow("Warning: Multiple 'Intercept' coefficients found, selecting the first...\n"))
      int_id <- int_id[1]
    }
    # If no explicit intercept, set as 1st value...
    if (length(int_id) == 0){
      int_id <- 1
    }
    base_val <- mod_obj$model$summary.fixed$mean[int_id]
    # if backtrans - backrans intercept to get baseline val
    if (!is.null(backtransform)){
      base_val <- backtransform(base_val)
    }
    
    
    # Get some basic info on model covars
    # all_covars <- all.vars(formula(delete.response(terms(obj$inla_fitted_formula))))
    all_covars <- all.vars(formula(delete.response(terms(mod_obj$user_formula$full))))
    # can ignore iid ranefs...
    # but covars are all needed... - can set to NA if missing though
    
    # TODO: put into hib_mod_proj_inla
    
    # Check type of data 
    message(crayon::green("Preparing input data for projection...\n"))
    if (is.data.frame(data)){
      # Don't need to do anything...Just copy data to data_df...
      data_df <- data
    }
    
    #   get from terr_proj_prep...
    # need to know fix_vars for terra_proj style approach...
    fix_vars <- all.vars(formula(delete.response(terms(mod_obj$user_formula$response_fixed))))
    # ignore "Intercept" too
    fix_vars <- fix_vars[!fix_vars == "Intercept"]
    
    # Get random vars
    rand_vars <- all.vars(mod_obj$user_formula$random)
    
    # And re_form_vars
    re_form_vars <- NA
    if (inherits(inla_specs$re_form, "formula")){
      # print("Checking random variables")
      re_form_vars <- all.vars(inla_specs$re_form)
      assertthat::assert_that(all(re_form_vars %in% rand_vars),
                              msg = "Specified re_form structure contains variables not in the original model.")
    }
    if (is.null(inla_specs$re_form)){
      re_form_vars <- rand_vars
    }

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
        # TODO: is ".," needed?
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
      
      # A <- inla.spde.make.A(input_mesh, 
      #                       loc = coords)
      s.site <- INLA::inla.spde.make.index('s', n.spde = site.spde$n.spde)
      # Create A, X, stk for prediction of dropped pop - using same mesh and spde?
      if (("Longitude" %in% colnames(data_df)) &
          ("Latitude" %in% colnames(data_df)) ){
        coords <- as.matrix(data_df[,c("Longitude", "Latitude")])
      }
      if (("Longitude_midpoint" %in% colnames(data_df)) &
          ("Latitude_midpoint" %in% colnames(data_df)) ){
        coords <- as.matrix(data_df[,c("Longitude_midpoint", "Latitude_midpoint")])
      }
      if (("x" %in% colnames(data_df)) &
          ("y" %in% colnames(data_df)) ){
        coords <- as.matrix(data_df[,c("x", "y")])
      }
      if (("X" %in% colnames(data_df)) &
          ("Y" %in% colnames(data_df)) ){
        coords <- as.matrix(data_df[,c("X", "Y")])
      }
      
      if (mod_obj$inla_mesh_settings$mesh_dim == "3d"){
        coords <- INLA::inla.mesh.map(coords,  
                                      projection = "longlat")
      }
      # if (obj$inla_mesh_settings$mesh_dim == "2d"){
      #   
      # }
      
      pred_A <- INLA::inla.spde.make.A(mod_obj$inla_mesh,
                                       loc = coords)
      
    } 
    
    # # Set up X data for modelling
    # pred_X <- All cols in data-df, apart form lat/long
    # if Intercept missing, add
    # if any other cols missing, e.g. SS/SSB, add as NA
    covars_sub <- all_covars[all_covars %in% colnames(data_df)]
    

    # print(all_covars)
    # print(covars_sub)
    
    pred_X <- data_df |>
      dplyr::select(covars_sub)
    
    if ( ("Intercept" %in% all_covars) &
         (!"Intercept" %in% covars_sub)){
      pred_X$Intercept <- 1
    }
    
    missing_covars <- all_covars[!all_covars %in% colnames(pred_X)]
    for (covar in missing_covars){
      pred_X[, covar] <- NA
    }
    
    # TODO: for categorical vars, ensure levels in projection input are within 
    # model fitting data (or NA)
    # for each fix_vars, check if chr or fct
    # if it is, check all data[,var] are in obj$data[,var
    check_covar_levels(mod_obj$model, pred_X, fix_vars)
    
    # TODO: also need to check any random vars - difference between fix_vars and all_vars
    if (inherits(inla_specs$re_form, "formula")){
        # print("Checking random variables")
        
        # Need to set rand_vars that are non-re_form_vars to NA
        non_reform_rand <- rand_vars[!rand_vars %in% re_form_vars]
        pred_X[, non_reform_rand] <- NA
        check_covar_levels(mod_obj$model, pred_X, re_form_vars)
    }
    if (is.null(inla_specs$re_form)){
        # print("Checking random variables")
        check_covar_levels(mod_obj$model, pred_X, re_form_vars)
    }
    if (is.na(re_form_vars)){
       # Need to set rand_vars to NA
       pred_X[, rand_vars] <- NA
    }
    # print(re_form_vars)
    # print(head(pred_X))

    # print(missing_covars)
    # print(head(pred_X))
    # Combine data etc into stack
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
    comb_stk <- INLA::inla.stack(mod_obj$inla_stack,
                                 pred_stk)
    message(crayon::green("Re-fitting INLA model to make projection (this may take a while)...\n"))
    # re-fit/ Project
    # TODO: incl while to ensure fixef are fitted properly!
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
                   control.predictor = list(A = INLA::inla.stack.A(comb_stk)),
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
                 control.predictor = list(A = INLA::inla.stack.A(comb_stk)),
                 verbose = F)
    }
    # %%%%%

    # m0 <- INLA::inla(formula = mod_obj$inla_fitted_formula,
    #                  family = mod_obj$model$.args$family,
    #                  data= INLA::inla.stack.data(comb_stk),
    #                  control.fixed = list(expand.factor.strategy='inla'),
    #                  # control.inla = list(strategy = "", int.strategy = ""), 
    #                  control.compute = list(dic = TRUE, waic = TRUE),
    #                  control.predictor = list(A = INLA::inla.stack.A(comb_stk)),
    #                  verbose = F)
    
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
      pred_df <- backtransform(pred_df)
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
    if (rescale_pred){
      # divide predictions by baseline val
      pred_cols <- colnames(pred_df)[grepl("prediction_", colnames(pred_df))]
      pred_rescaled_df <- pred_df[,pred_cols] / base_val
      colnames(pred_rescaled_df) <- paste("rescaled", colnames(pred_rescaled_df), sep = "_")
      # pred_rescaled_df <- pred_rescaled_df |>
      #   dplyr::select(-rescaled_prediction_sd)
    }
    
    pred_df <- dplyr::bind_cols(pred_df, pred_rescaled_df)
    
    
    # make back into data... - if df
    if (is.data.frame(data)){
      data <- pred_df
    }
    pred_cols <- colnames(pred_df)[grepl("prediction_", colnames(pred_df))]
    
    # if input was terra/list of terra, convert pred df back to terra format
    if (inherits(data, "SpatRaster")){
      # TODO: crs of input...
      # sf::st_crs(data)
      # data_pred <- terra::rast(pred_df[,c("x", "y", pred_cols)])
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
                           # TODO: crs of input...
                           # sf::st_crs(data_orig)
                           tmp_r <- terra::rast(tmp_df)
                           tmp_r <- terra::rast(tmp_df,
                                                crs = terra::crs(data_orig))
                           
                           out <- c(data_orig[[.x]], tmp_r)
                         })
      
      names(data) <- as.character(unique(pred_df$Year))
    }
    
    # return(data)
    
    # %%%%%
    # BRMS
    # %%%%%
  } else if (inherits(mod_obj, "hib_brms_model")){
    message(crayon::green("brms model detected! Projection in progress...\n"))
    
    # %%%%
    # Get incl_obs_uncert info
    # list(incl_obs_uncert = F)
    # print(brms_specs)
    if (!is.null(brms_specs)){
      if ("incl_obs_uncert" %in% names(brms_specs)){
        assertthat::assert_that(is.logical(brms_specs$incl_obs_uncert))
      }
      # if absent, set to F
      if (!"incl_obs_uncert" %in% names(brms_specs)){
        brms_specs$incl_obs_uncert <- FALSE
      }

      if ("re_form" %in% names(brms_specs)){
        assertthat::assert_that(
          is.na(c(brms_specs$re_form)) || is.null(brms_specs$re_form) || inherits(brms_specs$re_form, "formula"))
      }
      # if absent, set to F
      if (!"re_form" %in% names(brms_specs)){
        message(crayon::yellow("Warning: re_form missing, using default value of NA...\n"))
        brms_specs$re_form <- NA
      }
  
    }
    if (is.null(brms_specs)){
      message(crayon::yellow(paste0("Warning: brms_specs missing, using default values:\n",
                                  "\tincl_obs_uncert = FALSE,\n",
                                  "\tre_form = NA",
                                  " ...\n")))
      brms_specs <- list(incl_obs_uncert = FALSE,
                          re_form = NA)
    }
    # print(brms_specs)
    # %%%%
    
    # Get baseline value for optional rescaling
    if ("Intercept" %in% mod_obj$fixedef_coefficients$Coefficient_name){
      base_val <- mod_obj$fixedef_coefficients$Coefficient_value[mod_obj$fixedef_coefficients$Coefficient_name == "Intercept"]
    } else {
      base_val <- mod_obj$fixedef_coefficients$Coefficient_value[1]
    }
    
    # if backtrans - backrans intercept to get baseline val
    if (!is.null(backtransform)){
      base_val <- backtransform(base_val)
    }
    
    message(crayon::green("Preparing input data for projection...\n"))
    # Get some basic info on model covars
    all_covars <- all.vars(formula(delete.response(terms(mod_obj$model$formula$formula))))
    # all_covars <- all.vars(formula(obj$model)[[3]])
    # print(all_covars)
    # But what about random effects... - 
    rand_vars <- names(brms::ranef(mod_obj$model))
    fix_vars <- all_covars[!all_covars %in% rand_vars]
    # print(fix_vars)
    
    # Project for possible data types
    # Dataframe 
    if (is.data.frame(data)){
      # Check all fixef vars are in colnames of data
      assertthat::assert_that(all(fix_vars %in% colnames(data)),
                              msg = "Fixed effect variables missing from input data.")
      
      
      # TODO: for categorical vars, ensure levels in projection input are within 
      # model fitting data
      # for each fix_vars, check if chr or fct
      # if it is, check all data[,var] are in obj$data[,var
      check_covar_levels(mod_obj$model, data, fix_vars)
      
      # TODO: also need to check any random vars in brms_specs$re_form
      # Check they are present, and check they have the correct levels...
      if (inherits(brms_specs$re_form, "formula")){
          # print("Checking random variables")
          re_form_vars <- all.vars(brms_specs$re_form)
          assertthat::assert_that(all(re_form_vars %in% rand_vars), 
                                  msg = "Specified re_form structure contains variables not in the original model.")
          # all.vars(brms_specs$re_form[[2]])
          check_covar_levels(mod_obj$model, data, re_form_vars)
      }
      if (is.null(brms_specs$re_form)){
        check_covar_levels(mod_obj$model, data, rand_vars)
      }

      # %%%%%
      if (brms_specs$incl_obs_uncert == TRUE){
      # if (incl_obs_uncert == T){
        mean_df <- predict(mod_obj$model,
                           data,
                           # re_form = NA
                           re_form = brms_specs$re_form) |>
          as.data.frame()
        med_df <- predict(mod_obj$model,
                          data,
                          # re_form = NA,
                          re_form = brms_specs$re_form,
                          robust = TRUE) |>
          as.data.frame()
      }
      if (brms_specs$incl_obs_uncert == FALSE){
      # if (incl_obs_uncert == F){
        # brms::fitted.brmsfit
        mean_df <- fitted(mod_obj$model,
                          data,
                          # re_form = NA
                          scale = "response", # brms_specs$scale,
                          re_form = brms_specs$re_form) |>
        # re_form = brms_specs$re_form
          as.data.frame()
          # brms::fitted.brmsfit
        med_df <- fitted(mod_obj$model,
                         data,
                         # re_form = NA,
                         scale = "response", # brms_specs$scale,
                         re_form = brms_specs$re_form,
                         robust = TRUE) |>
          as.data.frame()
      }
      mean_df <- mean_df |>
        dplyr::rename(prediction_mean  = Estimate,
                      prediction_sd    = Est.Error,
                      prediction_q2.5  = Q2.5,
                      prediction_q97.5 = Q97.5)
      med_df <- med_df |>
        dplyr::rename(prediction_median= Estimate,
                      prediction_mad   = Est.Error)
      
      pred_df <- mean_df |>
        dplyr::bind_cols(med_df |>
                           dplyr::select(prediction_median,
                                         prediction_mad))
      # %%%%%
      
      # Use predict
      # pred_df <- as.data.frame(predict(obj$model,
      #                                  data,
      #                                  re_form = NA))
      
      # fix names...
      # pred_df <- pred_df |>
      #   dplyr::rename(prediction       = Estimate,
      #                 prediction_error = Est.Error,
      #                 prediction_q2.5  = Q2.5,
      #                 prediction_q97.5 = Q97.5)
      
      # What if we want to include  random effect??
      
      if (!is.null(backtransform)){
        pred_df <- backtransform(pred_df)
      }
      # %%%
      # Optional rescaling
      pred_rescaled_df <- NULL
      if (rescale_pred){
        # divide predictions by baseline val
        pred_cols <- colnames(pred_df)[grepl("prediction_", colnames(pred_df))]
        pred_rescaled_df <- pred_df[,pred_cols] / base_val
        colnames(pred_rescaled_df) <- paste("rescaled", colnames(pred_rescaled_df), sep = "_")
      }
      # pred_rescaled_df <- pred_rescaled_df |>
      #   dplyr::select(-rescaled_prediction_error)
      
      data <- dplyr::bind_cols(data, pred_df, pred_rescaled_df)
    }
    
    # Single raster object
    if (inherits(data, "SpatRaster")){
      # print(brms_specs)
      data <- terra_proj(data, mod_obj, backtransform, rescale_pred, base_val, fix_vars,
                          rand_vars,
                         brms_specs = brms_specs)
                         # incl_obs_uncert)
    }
    # List of raster objects
    if (( is.list(data) & !is.data.frame(data) & all(unlist(lapply(data, function(x){inherits(x, "SpatRaster")}))) ) ){
      # & inherits(data[[1]], "SpatRaster")
      # Apply individual terra proj to each item of list...
      data <- purrr::map(data,
                         terra_proj,
                         obj = mod_obj,
                         backtransform = backtransform,
                         rescale_pred = rescale_pred,
                         base_val = base_val,
                         fix_vars = fix_vars,
                         rand_vars = rand_vars,
                         brms_specs = brms_specs)
                         # incl_obs_uncert = incl_obs_uncert)
    }
    
  }
  
  # create projection object - need model obj, projection output, and projection_specs (backrans, rescale)
  # message(crayon::green"Creating projection object...\n"))
  proj <- hib_mod_projection$new(obj = mod_obj,
                                 proj = data,
                                 backtransform = backtransform,
                                 rescale_pred = rescale_pred,
                                 lme4_specs = lme4_specs,
                                 inla_specs = inla_specs,
                                 brms_specs = brms_specs,
                                 aoa_res    = aoa_res)
  return(proj)
}


#' # A helper function to wrap the preparation of SpatRaster data, and then project 
#' a hib, mixed-effects model - for lme4 and brms models
#' 
#' @description 
#' This function projects a fitted hib model using the supplied model object and
#' terra data to be used in projection
#' @details 
#' Expects a [`SpatRaster`].
#' @param data A [`SpatRaster`] that can be used to project a model.
#' @param obj A [`hib_model`] object - either lme4 or inla.
#' @param backtransform A [`function`] that will backtransform the model response
#' to the scale of the original variable, e.g., exp() as the inverse of log()
#' @param rescale_pred A [`logical`] value - should the predictions be rescaled relative
#' to the model intercept (assumes intercept is for baseline environment/LU)
#' @param base_val Model intercept, for use in rescaling predictions
#' @param fix_vars Model fixed effects
#' @param rand_vars Model random effects
#' @param lme4_specs A list containing details relevant to projecting an lme4-style
#' model. Specifically, named items should consist of 'bootstrap' [`logical`], and 
#' 'n_boot' [`numeric`], indicating whether bootstrapped prediction intervals should be 
#' estimated, and if so, how many bootstrap samples should be used.
#' @param brms_specs A list containing details relevant to projecting a brms-style
#' model. Specifically, named item should consist of 'incl_obs_uncert' [`logical`] 
#' to determine whether or not observation-level uncertainty is included in predictions. 
#' Setting this parameter to TRUE will produce predictions with larger uncertainty.
#' @keywords model, projections
#' @examples 
#' \dontrun{
#'   tmp_proj <- terra_proj(obj, data, backtransform, rescale_pred, base_val, fix_vars)
#' } 
#' @noRd
terra_proj <- function(data, obj, backtransform, rescale_pred, base_val, fix_vars,
                        rand_vars,
                       lme4_specs = NULL,
                       brms_specs = NULL){
                       # incl_obs_uncert = F,
                       # lme4_boot = T,
                       # n_boot = 100){
  # TODO: maybe change to @NoRd
  # assertions
  assertthat::assert_that(
    inherits(obj, "hib_lme4_model") || inherits(obj, "hib_brms_model"),
    inherits(data, "SpatRaster"), 
    is.function(backtransform) || is.null(backtransform),
    is.logical(rescale_pred),
    is.numeric(base_val),
    is.character(fix_vars),
    is.character(rand_vars) || is.na(c(rand_vars)) || is.null(rand_vars),
    # is.logical(incl_obs_uncert),
    # is.logical(lme4_boot),
    # is.numeric(n_boot) || is.null(n_boot)
    is.list(lme4_specs) || is.null(lme4_specs),
    is.list(brms_specs) || is.null(brms_specs)
  )
  # Covar checks
  # print(fix_vars)
  # print(names(data))

  # Extract important re_form_vars
  re_form_vars <- NA
  if (inherits(obj, "hib_lme4_model")){
    if (inherits(lme4_specs$re_form, "formula")){
      re_form_vars <- all.vars(lme4_specs$re_form)
      assertthat::assert_that(all(re_form_vars %in% rand_vars), 
                              msg = "Specified re_form structure contains variables not in the original model.")
    }
    if (is.null(lme4_specs$re_form)){
      re_form_vars <- rand_vars
    }
  }
  if (inherits(obj, "hib_brms_model")){
    if (inherits(brms_specs$re_form, "formula")){
      re_form_vars <- all.vars(brms_specs$re_form)
      assertthat::assert_that(all(re_form_vars %in% rand_vars), 
                              msg = "Specified re_form structure contains variables not in the original model.")
    }
    if (is.null(brms_specs$re_form)){
      re_form_vars <- rand_vars
    }
  }


  if (all(fix_vars %in% names(data))){
    # data_pred <- terra::predict(data,
    #                             obj$model,
    #                             re.form = NA)
    # 
    # if (!is.null(backtransform)){
    #   data_pred <- backtransform(data_pred)
    # }
    # names(data_pred) <- "prediction"
    
    # TODO: for categorical vars, ensure levels in projection input are within 
    # model fitting data
    # for each fix_vars, check if chr or fct
    # if it is, check all data[,var] are in obj$data[,var
    check_covar_levels(obj$model, data, fix_vars)

    if (inherits(obj, "hib_lme4_model")){
      # TODO: Random covar checks
      if (inherits(lme4_specs$re_form, "formula")){
        # print("Checking random variables")
        # re_form_vars <- all.vars(lme4_specs$re_form)
        # all.vars(brms_specs$re_form[[2]])
        check_covar_levels(obj$model, data, re_form_vars)
      }
      if (is.null(lme4_specs$re_form)){
        check_covar_levels(obj$model, data, re_form_vars)
      }

      # TODO: incorp uncertainty
      data_pred <- terra_proj_lme4_pred(obj$model, data, backtransform, lme4_specs)#lme4_boot, n_boot)
    }
    if (inherits(obj, "hib_brms_model")){
      # TODO: change to fitted...??
      # TODO: Random covar checks
      if (inherits(brms_specs$re_form, "formula")){
        # print("Checking random variables")
        # re_form_vars <- all.vars(brms_specs$re_form)
        # all.vars(brms_specs$re_form[[2]])
        check_covar_levels(obj$model, data, re_form_vars)
      }
      if (is.null(brms_specs$re_form)){
        check_covar_levels(obj$model, data, re_form_vars)
      }
      # print(brms_specs)
      data_pred <- terra_proj_brms_pred(obj$model, data, backtransform, brms_specs)#incl_obs_uncert)
    }
    # %%%%
  } else {
    # Any layer with a name in fix_vars is good to use 
    covars_in_rast <- names(data)[names(data) %in% fix_vars]
    # TODO: want to incl ranefs here -- either from re_form, or model
    if (!is.na(re_form_vars)){
      ranefs_in_rast <- names(data)[names(data) %in% re_form_vars]
      covars_in_rast <- c(covars_in_rast, ranefs_in_rast)
      # print(covars_in_rast)
    }

    tmp_pred_rast <- NULL
    n_pred_rast_lyrs <- 0
    if (length(covars_in_rast) > 0 ){
      tmp_pred_rast <- data[[covars_in_rast]]  
      n_pred_rast_lyrs <- dim(tmp_pred_rast)[3]
    }
    
    # any layer with a name not in fix_covars needs to be aligned to a var
    non_covars_in_rast <- names(data)[!names(data) %in% fix_vars]
    # TODO: want to drop ranefs here -- either from re_form, or model ??
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
      # TODO: Should use model-fitting data from obj$model...
      # lmer: @frame
      # inla: $.args$data
      # brms: @frame
      
      # # lme4 fit data,
      if ((class(obj$model) == "lmerMod" || class(obj$model) == "glmerMod")){
        fit_dat <- obj$model@frame
      }
      # brms fit data
      if (class(obj$model) == "brmsfit"){
        fit_dat <- obj$model$data
      }
      # inla fit data
      if (class(obj$model) == "inla"){
        fit_dat <- obj$model$.args$data
      }
      # 
      covars_not_in_rast <- purrr::map_dfr(covars_not_in_rast,
                                           function(x){
                                             if (is.character(fit_dat[[x]])){
                                               # if (is.factor(obj$data[[x]])){
                                               tmp_out <- data.frame(covar_nm = x,
                                                                     covar_levs = unique(fit_dat[[x]]))
                                               # Could add covar_lev_idx here... - save doing later...
                                             } else if (is.factor(fit_dat[[x]])){
                                             # if (is.factor(obj$data[[x]])){
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
        dplyr::left_join(is_binary_layer) |>
        dplyr::group_by(covar_nm) |>
        dplyr::mutate(is_binary_covar = all(is_binary_layer, na.rm = T))
      
      # print(covars_not_in_rast)
      
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
        message(crayon::yellow("Warning: SpatRaster with multiple non-binary layers detected! Projections are not yet supported...\n"))
        stop("In development")
      }
      
      
      # %%%%%%%%%%%%%%%%%%%%
      # TODO: - if there is a non_covar in rast that doesn't have a matching covar/lev
      # should throw an error?
      # %%%%%%%%%%%%%%%%%%%%
      # print(covars_not_in_rast)
      # print(non_covars_in_rast)
      # ID layers that are not within covariates required
      problem_layers <- non_covars_in_rast[!non_covars_in_rast %in% covars_not_in_rast$covar_levs_rast]
      # print(problem_layers)
      # Are there other situations where a problem layer could arise - or actualy not be a problem??
      # Note: extra layers in input that are not relevant to the model - eg a non-existent
      # covar level (Not Pasture) is simply ignored!
      if (length(problem_layers) > 0){
        # Maybe better to throw a warning/error??
        message(crayon::red(paste0("Error: covariate layer/s in the projection ",
                                   "data do not align with the model fitting data!\n\n",
                                   "Please fix this, see details below:\n")
                            )
                )
        
  
        message(paste("Problematic layers: ",
                      paste(problem_layers, collapse = " ; "),
                      "\n") 
                )
        
        stop("Unknown covariate layer, see above.", call. = F)
        
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
        # print(tmp_prop_rast)
        # Want to end up with N raster stacks - 1 per prop level
        n_prop_covar_lyrs <- dim(tmp_prop_rast)[3]
        
        # print(tmp_prop_rast)
      }
      # Check tmp_prop_rast is a terra SpatRaster
      # Could optionally exit here if prepping for inla...
      
    }
    # Everything to here is also useful for inla...
    
    # check N layers of each
    # if n_prop_covar_lyrs > 1  - list of rasters - each named based on covar level
    if (n_prop_covar_lyrs > 1){
      
      non_bin_sub <- non_bin_covars |>
        dplyr::filter(!is.na(covar_levs_rast))
      
      # print(tmp_prop_rast)
      # print(non_bin_covars)
      # print(non_bin_sub)
      non_bin_r <- tmp_rast[[non_bin_sub$covar_levs_rast]]
      
      tmp_pred_rast_ls <- vector("list", n_prop_covar_lyrs)
      names(tmp_pred_rast_ls) <- rep("", n_prop_covar_lyrs)
      for (i in 1:n_prop_covar_lyrs){
        # TODO: need to be careful combining these, if first item is null,
        # output becomes a list of 1...
        
        tmp_ <- c(tmp_pred_rast, tmp_bin_rast, tmp_prop_rast[[i]])
        # If tmp is a list, get [[1]]
        if (is.list(tmp_)){
          tmp_ <- tmp_[[1]]
        }
        # print(tmp_)

        assertthat::assert_that(
          inherits(tmp_, "SpatRaster")
        )
        # If is a SpatRaser, use as is.

        tmp_pred_rast_ls[[i]] <- tmp_ # c(tmp_pred_rast, tmp_bin_rast, tmp_prop_rast[[i]])[[1]]
        
        # Check all covars in new SpatRaster
        # print(tmp_pred_rast_ls[[i]])

        # print(c(tmp_pred_rast, tmp_bin_rast, tmp_prop_rast[[i]]))
        # print(fix_vars)
        # print(names(tmp_pred_rast_ls[[i]]))

        assertthat::assert_that(all(fix_vars %in% names(tmp_pred_rast_ls[[i]])))
        
        # tmp_pred_rast_ls[[i]] <- terra_proj_lme4_pred(obj$model, 
        #                                               tmp_pred_rast_ls[[i]], 
        #                                               backtransform)
        
        # TODO: for categorical vars, ensure levels in projection input are within 
        # model fitting data
        # for each fix_vars, check if chr or fct
        # if it is, check all data[,var] are in obj$data[,var
        check_covar_levels(obj$model, tmp_pred_rast_ls[[i]], fix_vars)
        
        if (inherits(obj, "hib_lme4_model")){
          # TODO : ranef covar check
          if (inherits(lme4_specs$re_form, "formula")){
            check_covar_levels(obj$model, tmp_pred_rast_ls[[i]], re_form_vars)
          }
          if (is.null(lme4_specs$re_form)){
            check_covar_levels(obj$model, tmp_pred_rast_ls[[i]], re_form_vars)
          }
          
          tmp_pred_rast_ls[[i]] <- terra_proj_lme4_pred(obj$model, 
                                                        tmp_pred_rast_ls[[i]], 
                                                        backtransform,
                                                        lme4_specs)
        }
        if (inherits(obj, "hib_brms_model")){
          # TODO : ranef covar check
          if (inherits(brms_specs$re_form, "formula")){
            check_covar_levels(obj$model, tmp_pred_rast_ls[[i]], re_form_vars)
          }
          if (is.null(brms_specs$re_form)){
            check_covar_levels(obj$model, tmp_pred_rast_ls[[i]], re_form_vars)
          }
          tmp_pred_rast_ls[[i]] <- terra_proj_brms_pred(obj$model, 
                                                        tmp_pred_rast_ls[[i]], 
                                                        backtransform,
                                                        brms_specs)
        }
        
        
        # Get proportional weighting
        lyr_nm <- terra::cats(tmp_prop_rast)[[i]] |> 
          dplyr::filter(ID == as.vector(unname(unique(terra::values(tmp_prop_rast[[i]])))) ) |> 
          dplyr::pull(category)
        
        # print(tmp_prop_rast[[i]])
        # print(terra::cats(tmp_prop_rast)[[i]])
        # print(unique(terra::values(tmp_prop_rast[[i]])))
        
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
      
      # TODO: for categorical vars, ensure levels in projection input are within 
      # model fitting data
      # for each fix_vars, check if chr or fct
      # if it is, check all data[,var] are in obj$data[,var
      check_covar_levels(obj$model, tmp_pred_rast, fix_vars)
      
      # TODO : ranef covar check

      # Predict
      if (inherits(obj, "hib_lme4_model")){
        if (inherits(lme4_specs$re_form, "formula")){
            check_covar_levels(obj$model, tmp_pred_rast, re_form_vars)
        }
        if (is.null(lme4_specs$re_form)){
          check_covar_levels(obj$model, tmp_pred_rast, re_form_vars)
        }

        data_pred <- terra_proj_lme4_pred(obj$model, tmp_pred_rast, backtransform,
                                          lme4_specs)
                                          # lme4_boot, n_boot)
      }
      if (inherits(obj, "hib_brms_model")){
        if (inherits(brms_specs$re_form, "formula")){
            check_covar_levels(obj$model, tmp_pred_rast, re_form_vars)
        }
        if (is.null(brms_specs$re_form)){
          check_covar_levels(obj$model, tmp_pred_rast, re_form_vars)
        }
        data_pred <- terra_proj_brms_pred(obj$model, tmp_pred_rast, backtransform,
                                          brms_specs)
                                          # incl_obs_uncert)
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
    # data_pred_rescaled <- data_pred_rescaled[[names(data_pred_rescaled) != "prediction_error"]]
    # data_pred_rescaled <- data_pred_rescaled[[names(data_pred_rescaled) != "prediction_Est.error"]]
    
    names(data_pred_rescaled) <- paste0("rescaled_", names(data_pred_rescaled))
  }
  
  data <- c(data, data_pred, data_pred_rescaled)
  # return SpatRaster or a list of SpatRasters
  return(data)
}


#' # A helper function to prepare SpatRaster data for projection
#' 
#' @description 
#' This function prepares the supplied SpatRaster to be ready for model projection
#' @details 
#' Expects a [`SpatRaster`].
#' @param covar_nm A [`character`], the name of the model covariate 
#' @param covar_df A [`data.frame`], with info on covariate, levels, and the associated
#' names present in the raster data
#' @param data_rast A [`SpatRaster`] that can be used to project a model.
#' @param mod_data A [`data.frame`], or list (inla), of the data used to fit the model
#' @param input_type A [`character`], either "binary" or "non-binary", indicating 
#' the type of data in the raster - i.e., is each grid cell only a single level of 
#' a covariate, or can it contain proportions of multiple levels
#' @keywords data, preparation, model, projection
#' @examples 
#' \dontrun{
#'   r <- terra_covar_prep(covar_nm, covar_df, data_rast, mod_data, input_type = "binary")
#' } 
#' @noRd
terra_covar_prep <- function(covar_nm, covar_df, data_rast, mod_data, input_type = "binary"){
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
      message(crayon::yellow("Warning: Binary SpatRaster contains cells with >1 category! Layers will be discarded ...\n"))
      return (NULL)
    }  
  }
  
  # get variable level IDs
  # What if mod_data[[covar_nm]] is not a factor!!?? - manually set... 
  # - or, should this be dealt with previously! - i.e in hib_fit... - for brms
  # unique(mod_data[[covar_nm]])
  # 1:dplyr::n()
  fac_lev_df <- data.frame(ID = as.numeric(mod_data[[covar_nm]]),
                           category = mod_data[[covar_nm]]) |>
    unique()
  # print(fac_lev_df)
  
  
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
    
    # TODO: is ".," needed in terra::categories ??
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
    
    # print(r |>
    #         max() |>
    #         terra::as.factor())
    # print(fac_lev_df)
    
    r <- r |>
      max() |>
      terra::as.factor() |>
      terra::categories(
                        # layer = 0,
                        value = fac_lev_df)
    
  }
  
  # %%%
  names(r) <- rep(covar_nm, dim(r)[3]) 
  
  # print(r)
  return(r)
}


#' A helper function to project a model using SpatRaster data
#' 
#' @description 
#' This function projects/predicts the supplied model using the SpatRaster 
#' @details 
#' Expects a hib lme4 model and a [`SpatRaster`].
#' @param mod A fitted hib_lme4_mod object 
#' @param data A [`SpatRaster`], prepared for projection
#' @param backtransform A [`function`] that can backtransform the response variable
#' to its original/natual scale.
#' @param lme4_specs A list containing details relevant to projecting an lme4-style
#' model. Specifically, named items should consist of 'bootstrap' [`logical`], and 
#' 'n_boot' [`numeric`], indicating whether bootstrapped prediction intervals should be 
#' estimated, and if so, how many bootstrap samples should be used.
#' Additionally, 're_form', as 're.form' in lme4:::predict.merMod(), used to specify random
#' effects to condition on when predicting an lme4 model (formula, NULL, or NA).
#' @keywords model, projection
#' @examples 
#' \dontrun{
#'   pred <- terra_proj_lme4_pred(mod, data, backtransform)
#' } 
#' @noRd
terra_proj_lme4_pred <- function(mod, data, backtransform, 
                                 lme4_specs){
                                 # lme4_boot, n_boot){
  # TODO: maybe change to @NoRd
  # Assertions
  assertthat::assert_that(
    inherits(mod, "merMod"),
    inherits(data, "SpatRaster"), 
    is.null(backtransform) || is.function(backtransform),
    # is.logical(lme4_boot), 
    # is.numeric(n_boot) || is.null(n_boot)
    is.list(lme4_specs) || is.null(lme4_specs)
  )
  # Check all fixef covars are in data - already done before func call!
  
  # %%%
  # Project
  # TODO: could also incl df.. - more general
  

  # %%%%%
  # Get info re bootstrapping
  # Should already be checked!!! - CHECK
  # Maybe just include for safety
  # if (!is.null(lme4_specs)){
  #   # Check for bootstrap - logical
  #   if ("bootstrap" %in% names(lme4_specs)){
  #     assertthat::assert_that(is.logical(lme4_specs$bootstrap))
  #   }
  #   # if absent, set to F
  #   if (!"bootstrap" %in% names(lme4_specs)){
  #     lme4_specs$bootstrap <- FALSE
  #   }
    
    
  #   # Check for n_boot - NA/integer/numeric
  #   if (lme4_specs$bootstrap){
  #     if ("n_boot" %in% names(lme4_specs)){
  #       assertthat::assert_that(is.numeric(lme4_specs$n_boot))
  #     }
  #     # if absent, set to F
  #     if (!"n_boot" %in% names(lme4_specs)){
  #       message(crayon::yellow("Warning: n_boot missing, using default value of 50...\n"))
  #       lme4_specs$n_boot <- 50
  #     }
  #   }
  #   # if absent, if bootstrap, set to 50, if not bootstrap, set to NA
  #   if (lme4_specs$bootstrap == F){
  #     lme4_specs$n_boot <- NA
  #   }
  # }
  # if (is.null(lme4_specs)){
  #   lme4_specs <- list(bootstrap = F, n_boot = NA) # , re_form = NA
  # }
  
  # %%%
  if (!is.null(lme4_specs)){
    # Check for bootstrap - logical
    if ("bootstrap" %in% names(lme4_specs)){
      assertthat::assert_that(is.logical(lme4_specs$bootstrap))
    }
    # if absent, set to F
    if (!"bootstrap" %in% names(lme4_specs)){
      lme4_specs$bootstrap <- FALSE
    }
    
    
    # Check for n_boot - NA/integer/numeric
    if (lme4_specs$bootstrap){
      if ("n_boot" %in% names(lme4_specs)){
        assertthat::assert_that(is.numeric(lme4_specs$n_boot))
        # Could fix n_boot to 50 if is.na(n_boot) but bootstrap == TRUE
      }
      # if absent, set to F
      if (!"n_boot" %in% names(lme4_specs)){
        message(crayon::yellow("Warning: n_boot missing, using default value of 50...\n"))
        lme4_specs$n_boot <- 50
      }
    }
    # if absent, if bootstrap, set to 50, if not bootstrap, set to NA
    if (lme4_specs$bootstrap == FALSE){
      lme4_specs$n_boot <- NA
    }

    if ("re_form" %in% names(lme4_specs)){
      assertthat::assert_that(
        is.na(c(lme4_specs$re_form)) || is.null(lme4_specs$re_form) || inherits(lme4_specs$re_form, "formula")
      )
    }
    if (!"re_form" %in% names(lme4_specs)){
      message(crayon::yellow("Warning: re_form missing, using default value of NA...\n"))
      lme4_specs$re_form <- NA
    }

  }
  if (is.null(lme4_specs)){
    message(crayon::yellow(paste0("Warning: lme4_specs missing, using default values:\n",
                                  "\tbootstrap = FALSE,\n",
                                  "\tn_boot = NA,\n",
                                  "\tre_form = NA",
                                  " ...\n")))
    lme4_specs <- list(bootstrap = FALSE, n_boot = NA, re_form = NA) #
  }

  data_pred <- terra::predict(data,
                              mod,
                              type = "response", # lme4_specs$type, ??
                              re.form = lme4_specs$re_form, 
                              na.rm = T)
  names(data_pred) <- "prediction_mean"
  
  # %%%
  # print(lme4_specs)
  # %%%%%
  if (lme4_specs$bootstrap == TRUE){
  # if (lme4_boot == T){
    # %%%%
    message(crayon::green(paste("\nGenerating bootstrapped prediciton intervals using",
                                lme4_specs$n_boot, "bootstrap samples...")))
    
    
    tmp_rast_df <- terra::as.data.frame(data, xy = T)
    # Progress bar!!
    # pb <- txtProgressBar(min = 1, max = n_boot, style = 3)
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
    
    prediction_sd <- apply(boot_obj$t, 2, sd)
    prediction_median <- apply(boot_obj$t, 2, function(x) as.numeric(quantile(x, probs = 0.5, na.rm = TRUE)))
    prediction_q2.5 <- apply(boot_obj$t, 2, function(x) as.numeric(quantile(x, probs = 0.025, na.rm = TRUE)))
    prediction_q97.5 <- apply(boot_obj$t, 2, function(x) as.numeric(quantile(x, probs= 0.975, na.rm = TRUE)))
    # TODO:
    # Above apply approach is v slow for large data sets...
    # 
    
    uncert_df <- data.frame(prediction_median= prediction_median,
                            prediction_sd    = prediction_sd,
                            prediction_q2.5  = prediction_q2.5,  
                            prediction_q97.5 = prediction_q97.5)
    uncert_df <- dplyr::bind_cols(tmp_rast_df[,c("x", "y")], uncert_df)
    # TODO: crs of input...
    # data_pred_uncert <- terra::rast(uncert_df)
    data_pred_uncert <- terra::rast(uncert_df,
                                    crs = terra::crs(data))
    
    data_pred <- c(data_pred, data_pred_uncert)
    # %%%%
  }
  
  
  # %%%
  if (!is.null(backtransform)){
    data_pred <- backtransform(data_pred)
  }
  # %%% 
  return(data_pred)
}



# %%%%
#' A helper function to project a model using SpatRaster data
#' 
#' @description 
#' This function projects/predicts the supplied model using the SpatRaster 
#' @details 
#' Expects a hib lme4 model and a [`SpatRaster`].
#' @param mod A fitted hib_brms_mod object 
#' @param data A [`SpatRaster`], prepared for projection
#' @param backtransform A [`function`] that can backtransform the response variable
#' to its original/natual scale.
#' @param brms_specs A list containing details relevant to projecting a brms-style
#' model. Specifically, named item should consist of 'incl_obs_uncert' [`logical`] 
#' to determine whether or not observation-level uncertainty is included in predictions.
#' Setting this parameter to TRUE will produce predictions with larger uncertainty.
#' @keywords model, projection
#' @examples 
#' \dontrun{
#'   pred <- terra_proj_brms_pred(mod, data, backtransform)
#' } 
#' @noRd
terra_proj_brms_pred <- function(mod, data, backtransform, brms_specs){#incl_obs_uncert){
  # TODO: maybe change to @NoRd
  
  # Assertions
  assertthat::assert_that(
    inherits(mod, "brmsfit"),
    inherits(data, "SpatRaster"), 
    is.null(backtransform) || is.function(backtransform),
    # is.logical(incl_obs_uncert)
    is.list(brms_specs) || is.null(brms_specs)
  )
  # Check type of inputs
  # Check all fixef covars are in data 
  
  # TODO: for categorical vars, ensure levels in projection input are within 
  # model fitting data
  # for each fix_vars, check if chr or fct
  # if it is, check all data[,var] are in obj$data[,var
  
  # %%%%
  # Get incl_obs_uncert info
  # list(incl_obs_uncert = F)
  # print(brms_specs)
  
  if (!is.null(brms_specs)){
    if ("incl_obs_uncert" %in% names(brms_specs)){
      assertthat::assert_that(is.logical(brms_specs$incl_obs_uncert))
    }
    # if absent, set to F
    if (!"incl_obs_uncert" %in% names(brms_specs)){
      brms_specs$incl_obs_uncert <- FALSE
    }
    
    if ("re_form" %in% names(brms_specs)){
      assertthat::assert_that(
        is.na(c(brms_specs$re_form)) || is.null(brms_specs$re_form) || inherits(brms_specs$re_form, "formula"))
    }
    # if absent, set to F
    if (!"re_form" %in% names(brms_specs)){
      message(crayon::yellow("Warning: re_form missing, using default value of NA...\n"))
      brms_specs$re_form <- NA
    }

  }
  if (is.null(brms_specs)){
    message(crayon::yellow(paste0("Warning: brms_specs missing, using default values:\n",
                                "\tincl_obs_uncert = FALSE,\n",
                                "\tre_form = NA",
                                " ...\n")))
    brms_specs <- list(incl_obs_uncert = FALSE,
                        re_form = NA)
  }
  # %%%%
  
  # %%%
  # Project
  # TODO: could also incl df.. - more general
  
  # print("HERE")
  # print(brms_specs$incl_obs_uncert)
  # include median for prediction (incl_obs_uncert)
  # %%%%%
  # %%%%%
  if (brms_specs$incl_obs_uncert == T){
  # if (incl_obs_uncert == T){
    mean_pred <- terra::predict(data,
                                mod,
                                # re_form = NA,
                                re_form = brms_specs$re_form,
                                na.rm = T)
    med_pred <- terra::predict(data,
                               mod,
                               # re_form = NA,
                               re_form = brms_specs$re_form,
                               robust = T,
                               na.rm = T)
    
    names(mean_pred)[names(mean_pred) == "Estimate"] <- "prediction_mean"
    names(mean_pred)[names(mean_pred) == "Est.Error"] <- "prediction_sd"
    names(mean_pred)[names(mean_pred) == "Q2.5"] <- "prediction_q2.5"
    names(mean_pred)[names(mean_pred) == "Q97.5"] <- "prediction_q97.5"
    
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
  # if (incl_obs_uncert == F){
    tmp_df <- terra::as.data.frame(data, xy = T)
    # brms::fitted.brmsfit
    mean_pred <- fitted(mod, tmp_df, 
                        # re_form = NA
                        scale = "response", # brms_specs$scale
                        re_form = brms_specs$re_form) |>
      as.data.frame()
    # brms::fitted.brmsfit
    median_pred <- fitted(mod, tmp_df, 
                          # re_form = NA,
                          scale = "response", # brms_specs$scale
                          re_form = brms_specs$re_form,
                          robust = T) |>
      as.data.frame()
    
    
    mean_pred <- mean_pred |>
      dplyr::rename(prediction_mean  = Estimate,
                    prediction_sd = Est.Error,
                    prediction_q2.5  = Q2.5,
                    prediction_q97.5 = Q97.5)
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
      data_pred <- backtransform(data_pred)
    }
  }
  
  # %%%%%
  # %%%%%
  # data_pred <- terra::predict(data,
  #                             mod,
  #                             re_form = NA)
  # if (!is.null(backtransform)){
  #   data_pred <- backtransform(data_pred)
  # }
  # 
  # names(data_pred)[names(data_pred) == "Estimate"] <- "prediction"
  # names(data_pred)[names(data_pred) == "Est.Error"] <- "prediction_error"
  # names(data_pred)[names(data_pred) == "Q2.5"] <- "prediction_q2.5"
  # names(data_pred)[names(data_pred) == "Q97.5"] <- "prediction_q97.5"
  # %%%%
  return(data_pred)
}


# terra_proj_inla_pred <- function(obj, data, backtransform, rescal_pred){
#   # Assertions...
#   
# }   

#' # A helper function to wrap the preparation of SpatRaster data and then return
#' a dataframe - for use in inla model projections
#' 
#' @description 
#' This function prepares a terra rast for use in INLA projections
#' Returns a dataframe
#' @details 
#' Expects a [`SpatRaster`].
#' @param data A [`SpatRaster`] that can be used to project a model.
#' @param obj A [`hib_model`] object -inla.
#' @param fix_vars Model fixed effects
#' @param rand_vars Model random effects
#' @param re_form_vars Model random effects to condition upon when predicting
#' @keywords model, projections
#' @examples 
#' \dontrun{
#'   tmp_df <- terra_to_df_inla(obj, data, fix_vars)
#' } 
#' @noRd
terra_to_df_inla <- function(data, obj, fix_vars, rand_vars, re_form_vars){
  # TODO: maybe change to @NoRd
  # TODO: maybe co-opt to use fitted for brms rather then predict - for terra input
  
  # assertions
  assertthat::assert_that(
    inherits(obj, "hib_inla_model"),
    inherits(data, "SpatRaster"), 
    is.character(fix_vars),
    is.character(rand_vars) || is.na(c(rand_vars)) || is.null(rand_vars),
    is.character(re_form_vars) || is.na(c(re_form_vars)) || is.null(re_form_vars)
  )
  # Covar checks
  
  # TODO: What about non-fix (random) variables

  if (all(fix_vars %in% names(data))){
    # TODO: also check for all ranefs... re_form_vars
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

          message(crayon::red(paste0("Error: random effect variable specified in ",
                                    "re_form is missing.\n")))
          message(crayon::red(paste("Missing covariates: ",
                                missing_covars_chr, sep = "")))
          stop()  
        }
        
      }
    
  } else {
    # Any layer with a name in fix_vars is good to use 
    covars_in_rast <- names(data)[names(data) %in% fix_vars]
    # TODO: also get ranefs. - assume not proportions/binary - as in model
    if (!is.na(re_form_vars)){
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
    # TODO: drop ranefs in data from non_covars_in_rast
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

      # TODO: Should use model-fitting data from obj$model...
      # lmer: @frame
      # inla: $.args$data
      # brms: @frame
      
      # # lme4 fit data,
      if ((class(obj$model) == "lmerMod" || class(obj$model) == "glmerMod")){
        fit_dat <- obj$model@frame
      }
      # brms fit data
      if (class(obj$model) == "brmsfit"){
        fit_dat <- obj$model$data
      }
      # inla fit data
      if (class(obj$model) == "inla"){
        fit_dat <- obj$model$.args$data
      }
      # if (is.factor(fit_dat[[x]]))
      covars_not_in_rast <- purrr::map_dfr(covars_not_in_rast,
                                           function(x){
                                             if (is.factor(fit_dat[[x]])){
                                             # if (is.factor(obj$data[[x]])){
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
        dplyr::left_join(is_binary_layer) |>
        dplyr::group_by(covar_nm) |>
        dplyr::mutate(is_binary_covar = all(is_binary_layer, na.rm = T))
      
      # Split binary and non binary layers
      bin_covars <- covars_not_in_rast |>
        dplyr::filter(is_binary_covar == T)
      
      non_bin_covars <- covars_not_in_rast |>
        dplyr::filter(is_binary_covar == F)
      
      if (length(unique(non_bin_covars$covar_nm)) > 1){
        message(crayon::yellow("Warning: SpatRaster with multiple non-binary layers detected! Projections are not yet supported...\n"))
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
        # print(tmp_prop_rast)
        # Want to end up with N raster stacks - 1 per prop level
        n_prop_covar_lyrs <- dim(tmp_prop_rast)[3]
      }
      # Check tmp_prop_rast is a terra SpatRaster
      # Could optionally exit here if prepping for inla...
      
    }
    # Everything to here is also useful for inla...
    
    # check N layers of each
    # if n_prop_covar_lyrs > 1  - list of rasters - each named based on covar level
    if (n_prop_covar_lyrs > 1){
      
      non_bin_sub <- non_bin_covars |>
        dplyr::filter(!is.na(covar_levs_rast))
      
      non_bin_r <- tmp_rast[[non_bin_sub$covar_levs_rast]]
      
      tmp_pred_df_ls <- vector("list", n_prop_covar_lyrs)
      names(tmp_pred_df_ls) <- rep("", n_prop_covar_lyrs)
      for (i in 1:n_prop_covar_lyrs){
        # TODO: need to be careful combining these, if first item is null,
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

        tmp_r <- tmp_ # c(tmp_pred_rast, tmp_bin_rast, tmp_prop_rast[[i]])[[1]]
        
        # %%%%
        # Check fix vars and re-form_vars are present
        assertthat::assert_that(all(fix_vars %in% names(tmp_r)))
        check_covar_levels(obj$model, tmp_r, fix_vars)
        
        if (!is.na(re_form_vars)){
          check_covar_levels(obj$model, tmp_r, re_form_vars)
        }        
        # %%%%

        # %%%
        # tmp_r <- c(tmp_pred_rast, tmp_bin_rast, tmp_prop_rast[[i]])[[1]]
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
        # print(head(prop_df))
        # rename lyr_orig_nm as prop_cover
        # left_join prop_df to tmp_df
        tmp_pred_df_ls[[i]] <- dplyr::left_join(tmp_df, prop_df[,c("x", "y", "prop_cover")])
        
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
  
  # print(head(data_df))
  return(data_df)  
  
}


#' An internal helper function to check that the data provided for the projection 
#' does not contain any factor-levels that are not in the model fitting data.
#' 
#' @description 
#' This function checks whether the supplied projection data is compatible with
#' the supplied model
#' Returns either nothing (NULL, input data passes checks), or throws an error 
#' (input data fails checks)
#' 
#' @details 
#' Expects a model object (lme4/brms/?inla?).
#' @param mod A fitted model, either from lme4, inla or brms
#' @param proj_dat The prepared data used to predict from
#' @param covars The fixed-effect covariates in the model # TODO: could also be ranefs
#' @keywords model projections data checks
#' @examples 
#' \dontrun{
#'   check_covar_levels(model, projection_df, fixed_covars)
#' } 
#' @noRd
check_covar_levels <- function(mod, proj_dat, covars){
  
  # Assertion checks
  assertthat::assert_that(
    class(mod) == "lmerMod" || class(mod) == "glmerMod" || 
      class(mod) == "brmsfit" || class(mod) == "inla",
    is.data.frame(proj_dat) || inherits(proj_dat, "SpatRaster"),
    is.character(covars)
  )
  
  # First check - that all covars are in proj_dat
  if (!all(covars %in% names(proj_dat))){
    
    # Get missing covars
    missing_covars <- covars[!covars %in% names(proj_dat)]
    missing_covars_chr <- paste(missing_covars, collapse = " ; ")
    message(crayon::red("Error: covariates missing from projection data!\n"))
    message(crayon::red(paste("Missing covariates: ",
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
  if ((class(mod) == "lmerMod" || class(mod) == "glmerMod")){
    fit_dat <- mod@frame
  }
  # brms fit data
  if (class(mod) == "brmsfit"){
    fit_dat <- mod$data
  }
  # inla fit data
  if (class(mod) == "inla"){
    fit_dat <- mod$.args$data
  }
  
  problem_vars <- list()
  
  # print(fit_dat)
  
  is_problem_var <- FALSE
  for (i in 1:length(covars)){
    # Check if var is a factor/chr
    var_i <- covars[i]
    # print(var_i)
    if ( (is.factor(fit_dat[[var_i]]) | is.character(fit_dat[[var_i]])) ){
      # Get unique chr values for fit data
      var_i_fit_levs <- fit_dat[[var_i]] |> 
        unique() |>
        as.character()
      # Alternative - get levels directly - if factor
      # levels(fit_dat[[var_i]])
      if (is.factor(fit_dat[[var_i]])){
        var_i_fit_levs <- levels(fit_dat[[var_i]])
      }
      
      # print(var_i_fit_levs)
      # Get unique chr values for proj data 
      var_i_proj <- proj_dat[[var_i]]
      
      if (proj_dat_type == "df"){
        var_i_proj_levs <- var_i_proj |> 
          unique() |>
          as.character()
        # Alternative - get levels directly - if factor
        # levels(fit_dat[[var_i]])
        if (is.factor(var_i_proj)){
          var_i_proj_levs <- levels(var_i_proj)
        }
      }
      if (proj_dat_type == "rast"){
        var_i_proj_levs <- terra::levels(var_i_proj)[[1]][["category"]]
      }
      
      # print(head(var_i_proj_levs))
      # print(head(var_i_fit_levs))
      # print(var_i_proj_levs %in% var_i_fit_levs)

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
    message(crayon::red(paste0("Some categorical covariates in the projection data ",
                               "are not in the model fitting data!\n\n",
                               "Please fix this, see details below:\n"
                               ) ))
    problem_vars |>
      dplyr::bind_rows() |>
      insight::format_table(use_symbols = T, digits = 0) |>
      insight::export_table() |>
      message()
    
    # print()
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






