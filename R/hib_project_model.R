#' Project a hib, mixed-effects model
#' 
#' @description 
#' This function projects a fitted hib model using the supplied model object and
#' data to be used in projection - i.e., a df or raster stack of covars
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`hib_model`] object - either lme4 or inla.
#' @param data Either a dataframe-like object or a raster stack that can be used
#' to project onto
#' @param backtransform A [`function`] that will backtransform the model response
#' to the scale of the original variable, e.g., exp() as the inverse of log()
#' @param rescale_pred A [`logical`] value - should the predictions be rescaled relative
#' to the model intercept (assumes intercept is for baseline environment/LU)
#' @param incl_obs_uncert A [`logical`] indicating whether or not observation-level
#' uncertainty should be included in the predictions from brms models. Setting 
#' this parameter to TRUE will produce predictions with larger uncertainty.
#' @keywords model, projections
#' @examples 
#' \dontrun{
#'   proj <- hib_mod_project(obj, data)
#' } 
#' @export
hib_mod_project <- function(obj, 
                            data,
                            backtransform = NULL,
                            rescale_pred = F,
                            incl_obs_uncert = F){
  # Assertions 
  # is data a df/tibble/ raster/stack 
  assertthat::assert_that(
    inherits(obj, "hib_lme4_model") || inherits(obj, "hib_inla_model") || inherits(obj, "hib_brms_model"), # obj is a hib model
    is.data.frame(data) || inherits(data, "SpatRaster") || ( is.list(data) & !is.data.frame(data) & all(unlist(lapply(data, function(x){inherits(x, "SpatRaster")}))) ), 
    is.function(backtransform) || is.null(backtransform),
    is.logical(rescale_pred),
    is.logical(incl_obs_uncert)
  )
  # TODO: better check of input data - all items in list ar SpatRasters
  # all(unlist(lapply(data, function(x){inherits(x, "SpatRaster")})))
  
  # Check type of model
  # %%%%%
  # LME4
  # %%%%%
  if (inherits(obj, "hib_lme4_model")){
    message(crayon::green("lme4 model detected! Projection in progress...\n"))
    # TODO: wrap in hib_mod_proj_lme4??
    
    # Get baseline value for optional rescaling
    
    # %%%%
    # get model intercept
    if ("(Intercept)" %in% obj$fixedef_coefficients$Coefficient_name){
      base_val <- obj$fixedef_coefficients$Coefficient_value[obj$fixedef_coefficients$Coefficient_name == "(Intercept)"]
    } else {
      base_val <- obj$fixedef_coefficients$Coefficient_value[1]
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
    all_covars <- all.vars(formula(delete.response(terms(formula(obj$model)))))
    # But what about random effects... - 
    rand_vars <- names(lme4::ranef(obj$model))
    fix_vars <- all_covars[!all_covars %in% rand_vars]
    
    # Project for possible data types
    # Dataframe 
    if (is.data.frame(data)){
      # Check all fixef vars are in colnames of data
      assertthat::assert_that(all(fix_vars %in% colnames(data)))
      
      message(crayon::green("Projecting model...\n"))
      pred_df <- as.data.frame(predict(obj$model,
                                       data,
                                       re.form = NA))
      colnames(pred_df) <- "prediction"
      
      # %%%
      # TODO: incorp uncertainty
      # Use bootMer - 
      # or sampling from varCovar...??
      message(crayon::green("Bootstrapping model predictions...\n"))
      boot <- lme4::bootMer(obj$model, 
                            function(mod, df = data){predict(mod, df, re.form = NA)}, 
                            nsim = 100)
      
      pred_df$prediction_median <- apply(boot$t, 2, function(x) as.numeric(quantile(x, probs = 0.5, na.rm = TRUE)))
      pred_df$prediction_sd <- apply(boot$t, 2, sd)
      pred_df$prediction_q2.5 <- apply(boot$t, 2, function(x) as.numeric(quantile(x, probs = 0.025, na.rm = TRUE)))
      pred_df$prediction_q97.5 <- apply(boot$t, 2, function(x) as.numeric(quantile(x, probs = 0.975, na.rm = TRUE)))
      # std.err <- apply(boot_test$t, 2, sd)
      # CI.lower1 <- pred - std.err*1.96
      # CI.upper1 <- pred + std.err*1.96
      # %%%
      
      # What if we want to include  random effect??
      
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
      data <- terra_proj(data, obj, backtransform, rescale_pred, base_val, fix_vars)
    }
    # List of raster objects
    if ((is.list(data) & !is.data.frame(data))){
      # & inherits(data[[1]], "SpatRaster")
      # Apply individual terra proj to each item of list...
      data <- purrr::map(data,
                         terra_proj,
                         obj = obj, 
                         backtransform = backtransform, 
                         rescale_pred = rescale_pred, 
                         base_val = base_val, 
                         fix_vars = fix_vars)
    }
    
  } 
  # %%%%%
  # INLA
  # %%%%%
  else if (inherits(obj, "hib_inla_model")){
    # message(crayon::yellow("Warning: INLA model detected! Projections are not yet supported...\n"))
    message(crayon::green("INLA model detected! Projection in progress...\n"))
    # stop("In development")
    
    # Get baseline value for optional rescaling
    # get model intercept
    int_id <- which(grepl("intercept", rownames(obj$model$summary.fixed), ignore.case = T))
    # TODO: What if another covar also contains intercept...
    if (length(int_id) > 1){
      message(crayon::yellow("Warning: Multiple 'Intercept' coefficients found, selecting the first...\n"))
      int_id <- int_id[1]
    }
    # If no explicit intercept, set as 1st value...
    if (length(int_id) == 0){
      int_id <- 1
    }
    base_val <- obj$model$summary.fixed$mean[int_id]
    # if backtrans - backrans intercept to get baseline val
    if (!is.null(backtransform)){
      base_val <- backtransform(base_val)
    }
    
    
    # Get some basic info on model covars
    # all_covars <- all.vars(formula(delete.response(terms(obj$inla_fitted_formula))))
    all_covars <- all.vars(formula(delete.response(terms(obj$user_formula$full))))
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
    fix_vars <- all.vars(formula(delete.response(terms(obj$user_formula$response_fixed))))
    # ignore "Intercept" too
    fix_vars <- fix_vars[!fix_vars == "Intercept"]
    
    # TODO: check and test
    # If terra, convert rasters to dfs
    if (inherits(data, "SpatRaster")){
      data_df <- terra_to_df_inla(data, obj, fix_vars)
    }
    # TODO
    # If list of terra, convert each raster to df, add name of item in list, bind all
    if ((is.list(data) & !is.data.frame(data) &
         all(unlist(lapply(data, function(x){inherits(x, "SpatRaster")}))) ) ){
      data_df <- purrr::map(data,
                            terra_to_df_inla,
                            obj = obj, 
                            fix_vars = fix_vars) %>%
        # TODO: is ".," needed?
        # TODO: better term for id
        dplyr::bind_rows(., .id = "Year")
    }
    
    # With data_df...
    # need to check presence of covars, add NA if any missing, add 1 for intercept
    # and lat/long if spatial ranef...
    # obj$inla_fitted_formula
    # all_covars
    # obj$inla_incl_spatial_term
    # obj$inla_stack
    
    
    # %%%%
    
    if (obj$inla_incl_spatial_term){
      # Need s.spde from original model!!
      # or just recreate!
      site.spde <- INLA::inla.spde2.matern(obj$inla_mesh, alpha = 2)
      
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
      
      if (obj$inla_mesh_settings$mesh_dim == "3d"){
        coords <- INLA::inla.mesh.map(coords,  
                                      projection = "longlat")
      }
      # if (obj$inla_mesh_settings$mesh_dim == "2d"){
      #   
      # }
      
      pred_A <- INLA::inla.spde.make.A(obj$inla_mesh,
                                       loc = coords)
      
    } 
    
    # # Set up X data for modelling
    # pred_X <- All cols in data-df, apart form lat/long
    # if Intercept missing, add
    # if any other cols missing, e.g. SS/SSB, add as NA
    covars_sub <- all_covars[all_covars %in% colnames(data_df)]
    
    # print(all_covars)
    # print(covars_sub)
    
    pred_X <- data_df %>%
      dplyr::select(covars_sub)
    
    if ( ("Intercept" %in% all_covars) &
         (!"Intercept" %in% covars_sub)){
      pred_X$Intercept <- 1
    }
    
    missing_covars <- all_covars[!all_covars %in% colnames(pred_X)]
    for (covar in missing_covars){
      pred_X[,covar] <- NA
    }
    # print(missing_covars)
    # print(head(pred_X))
    # Combine data etc into stack
    if (obj$inla_incl_spatial_term){
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
    comb_stk <- INLA::inla.stack(obj$inla_stack,
                                 pred_stk)
    message(crayon::green("Re-fitting INLA model to make projection (this may take a while)...\n"))
    # re-fit/ Project
    m0 <- INLA::inla(formula = obj$inla_fitted_formula,
                     family = obj$model$.args$family,
                     data= INLA::inla.stack.data(comb_stk),
                     control.fixed = list(expand.factor.strategy='inla'),
                     # control.inla = list(strategy = "", int.strategy = ""), 
                     control.compute = list(dic = TRUE, waic = TRUE),
                     control.predictor = list(A = INLA::inla.stack.A(comb_stk)),
                     verbose = F)
    
    # Get pred_stack index
    pred_idx <- INLA::inla.stack.index(comb_stk, tag = "pred")$data
    
    pred_df <- m0$summary.fitted.values[pred_idx, ]
    rownames(pred_df) <- NULL
    
    colnames(pred_df) <- paste("prediction", colnames(pred_df), sep = "_")
    
    # TODO: need to harmonise prediction colnames
    
    # prediction_mean prediction_sd prediction_q2.5
    # prediction_0.5quant prediction_0.975quant prediction_mode
    
    # prediction_mean prediction_sd prediction_0.025quant
    # prediction_median prediction_q97.5 prediction_mode
    pred_df <- pred_df %>%
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
          pred_df <- pred_df %>%
            dplyr::group_by(x, y, Year) %>%
            dplyr::summarise(dplyr::across(dplyr::contains("prediction"), ~ stats::weighted.mean(.x, w = prop_cover)))
        }
        if (!"Year" %in% colnames(data_df)){
          pred_df <- pred_df %>%
            dplyr::group_by(x, y)  %>%
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
      # pred_rescaled_df <- pred_rescaled_df %>%
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
                           tmp_df <- df %>%
                             dplyr::filter(Year == .x) %>%
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
  } else if (inherits(obj, "hib_brms_model")){
    message(crayon::green("brms model detected! Projection in progress...\n"))
    
    # Get baseline value for optional rescaling
    if ("Intercept" %in% obj$fixedef_coefficients$Coefficient_name){
      base_val <- obj$fixedef_coefficients$Coefficient_value[obj$fixedef_coefficients$Coefficient_name == "Intercept"]
    } else {
      base_val <- obj$fixedef_coefficients$Coefficient_value[1]
    }
    
    # if backtrans - backrans intercept to get baseline val
    if (!is.null(backtransform)){
      base_val <- backtransform(base_val)
    }
    
    message(crayon::green("Preparing input data for projection...\n"))
    # Get some basic info on model covars
    all_covars <- all.vars(formula(delete.response(terms(obj$model$formula$formula))))
    # all_covars <- all.vars(formula(obj$model)[[3]])
    # print(all_covars)
    # But what about random effects... - 
    rand_vars <- names(brms::ranef(obj$model))
    fix_vars <- all_covars[!all_covars %in% rand_vars]
    # print(fix_vars)
    
    # Project for possible data types
    # Dataframe 
    if (is.data.frame(data)){
      # Check all fixef vars are in colnames of data
      assertthat::assert_that(all(fix_vars %in% colnames(data)))
      
      
      # TODO: use either predict() or fitted() - 
      # TODO: get median (robust = T) and mean (robust = F)
      # %%%%%
      if (incl_obs_uncert == T){
        mean_df <- predict(obj$model,
                           data,
                           re_form = NA) %>%
          as.data.frame()
        med_df <- predict(obj$model,
                          data,
                          re_form = NA,
                          robust = T) %>%
          as.data.frame()
      }
      if (incl_obs_uncert == F){
        mean_df <- fitted(obj$model,
                          data,
                          re_form = NA) %>%
          as.data.frame()
        med_df <- fitted(obj$model,
                         data,
                         re_form = NA,
                         robust = T) %>%
          as.data.frame()
      }
      mean_df <- mean_df %>%
        dplyr::rename(prediction_mean  = Estimate,
                      prediction_sd    = Est.Error,
                      prediction_q2.5  = Q2.5,
                      prediction_q97.5 = Q97.5)
      med_df <- med_df %>%
        dplyr::rename(prediction_median= Estimate,
                      prediction_mad   = Est.Error)
      
      pred_df <- mean_df %>%
        dplyr::bind_cols(med_df %>%
                           dplyr::select(prediction_median,
                                         prediction_mad))
      # %%%%%
      
      # Use predict
      # pred_df <- as.data.frame(predict(obj$model,
      #                                  data,
      #                                  re_form = NA))
      
      # fix names...
      # pred_df <- pred_df %>%
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
      # pred_rescaled_df <- pred_rescaled_df %>%
      #   dplyr::select(-rescaled_prediction_error)
      
      data <- dplyr::bind_cols(data, pred_df, pred_rescaled_df)
    }
    
    # Single raster object
    if (inherits(data, "SpatRaster")){
      data <- terra_proj(data, obj, backtransform, rescale_pred, base_val, fix_vars,
                         incl_obs_uncert)
    }
    # List of raster objects
    if (( is.list(data) & !is.data.frame(data) & all(unlist(lapply(data, function(x){inherits(x, "SpatRaster")}))) ) ){
      # & inherits(data[[1]], "SpatRaster")
      # Apply individual terra proj to each item of list...
      data <- purrr::map(data,
                         terra_proj,
                         obj = obj,
                         backtransform = backtransform,
                         rescale_pred = rescale_pred,
                         base_val = base_val,
                         fix_vars = fix_vars,
                         incl_obs_uncert = incl_obs_uncert)
    }
    
  }
  
  # create projection object - need model obj, projection output, and projection_specs (backrans, rescale)
  # message(crayon::green"Creating projection object...\n"))
  proj <- hib_mod_projection$new(obj = obj,
                                 proj = data,
                                 backtransform = backtransform,
                                 rescale_pred = rescale_pred)
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
#' @param incl_obs_uncert A [`logical`] indicating whether or not observation-level
#' uncertainty should be included in the predictions from brms models. Setting 
#' this parameter to TRUE will produce predictions with larger uncertainty.
#' @keywords model, projections
#' @examples 
#' \dontrun{
#'   tmp_proj <- terra_proj(obj, data, backtransform, rescale_pred, base_val, fix_vars)
#' } 
#' @export
terra_proj <- function(data, obj, backtransform, rescale_pred, base_val, fix_vars,
                       incl_obs_uncert = F){
  # TODO: maybe change to @NoRd
  # assertions
  assertthat::assert_that(
    inherits(obj, "hib_lme4_model") || inherits(obj, "hib_brms_model"),
    inherits(data, "SpatRaster"), 
    is.function(backtransform) || is.null(backtransform),
    is.logical(rescale_pred),
    is.numeric(base_val),
    is.character(fix_vars),
    is.logical(incl_obs_uncert)
  )
  # Covar checks
  # print(fix_vars)
  # print(names(data))
  if (all(fix_vars %in% names(data))){
    # data_pred <- terra::predict(data,
    #                             obj$model,
    #                             re.form = NA)
    # 
    # if (!is.null(backtransform)){
    #   data_pred <- backtransform(data_pred)
    # }
    # names(data_pred) <- "prediction"
    if (inherits(obj, "hib_lme4_model")){
      # TODO: incorp uncertainty
      data_pred <- terra_proj_lme4_pred(obj$model, data, backtransform)
    }
    if (inherits(obj, "hib_brms_model")){
      # TODO: change to fitted...??
      data_pred <- terra_proj_brms_pred(obj$model, data, backtransform, incl_obs_uncert)
    }
    # %%%%
  } else {
    # Any layer with a name in fix_vars is good to use 
    covars_in_rast <- names(data)[names(data) %in% fix_vars]
    tmp_pred_rast <- NULL
    n_pred_rast_lyrs <- 0
    if (length(covars_in_rast) > 0 ){
      tmp_pred_rast <- data[[covars_in_rast]]  
      n_pred_rast_lyrs <- dim(tmp_pred_rast)[3]
    }
    
    # any layer with a name not in fix_covars needs to be aligned to a var
    non_covars_in_rast <- names(data)[!names(data) %in% fix_vars]
    tmp_rast <- NULL
    if (length(non_covars_in_rast) > 0 ){
      tmp_rast <- data[[non_covars_in_rast]]
    }
    
    # Assume that layer names match levels of covariates
    # e.g "LUIPrimary_minimal use" or just "Primary_minimal use"
    covars_not_in_rast <- fix_vars[!fix_vars %in% names(data)]
    
    if (length(non_covars_in_rast) > 0 ){
      # for each covar, get unique levels (if factor in mod_df)
      covars_not_in_rast <- purrr::map_dfr(covars_not_in_rast,
                                           function(x){
                                             if (is.factor(obj$data[[x]])){
                                               tmp_out <- data.frame(covar_nm = x,
                                                                     covar_levs = levels(obj$data[[x]]))
                                               # Could add covar_lev_idx here... - save doing later...
                                             }
                                             else {
                                               tmp_out <- data.frame(covar_nm = x,
                                                                     covar_levs = NA)
                                             }
                                             tmp_out
                                           })
      
      covars_not_in_rast <- covars_not_in_rast %>%
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
                                                     is_binary_layer = all(unique(terra::values(tmp_rast[[nm]])) %in% c(0,1))
                                          )
                                        })
      
      covars_not_in_rast <- covars_not_in_rast %>%
        dplyr::left_join(is_binary_layer) %>%
        dplyr::group_by(covar_nm) %>%
        dplyr::mutate(is_binary_covar = all(is_binary_layer, na.rm = T))
      
      # Split binary and non binary layers
      bin_covars <- covars_not_in_rast %>%
        dplyr::filter(is_binary_covar == T)
      
      non_bin_covars <- covars_not_in_rast %>%
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
                                   mod_data = obj$data,
                                   input_type = "binary") %>%
          unlist() %>%
          c() %>%
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
                                    mod_data = obj$data,
                                    input_type = "non-binary") %>%
          unlist() %>%
          c() %>%
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
      
      non_bin_sub <- non_bin_covars %>%
        dplyr::filter(!is.na(covar_levs_rast))
      
      non_bin_r <- tmp_rast[[non_bin_sub$covar_levs_rast]]
      
      tmp_pred_rast_ls <- vector("list", n_prop_covar_lyrs)
      names(tmp_pred_rast_ls) <- rep("", n_prop_covar_lyrs)
      for (i in 1:n_prop_covar_lyrs){
        # TODO: need to be careful combining these, if first item is null,
        # output becomes a list of 1...
        tmp_pred_rast_ls[[i]] <- c(tmp_pred_rast, tmp_bin_rast, tmp_prop_rast[[i]])[[1]]
        # Check all covars in new SpatRaster
        assertthat::assert_that(all(fix_vars %in% names(tmp_pred_rast_ls[[i]])))
        
        # tmp_pred_rast_ls[[i]] <- terra_proj_lme4_pred(obj$model, 
        #                                               tmp_pred_rast_ls[[i]], 
        #                                               backtransform)
        if (inherits(obj, "hib_lme4_model")){
          tmp_pred_rast_ls[[i]] <- terra_proj_lme4_pred(obj$model, 
                                                        tmp_pred_rast_ls[[i]], 
                                                        backtransform)
        }
        if (inherits(obj, "hib_brms_model")){
          tmp_pred_rast_ls[[i]] <- terra_proj_brms_pred(obj$model, 
                                                        tmp_pred_rast_ls[[i]], 
                                                        backtransform,
                                                        incl_obs_uncert)
        }
        
        
        # Get proportional weighting
        lyr_nm <- terra::cats(tmp_prop_rast)[[i]] %>% 
          dplyr::filter(ID == as.vector(unname(unique(terra::values(tmp_prop_rast[[i]])))) ) %>% 
          dplyr::pull(category)
        
        lyr_orig_nm <- non_bin_sub %>%
          dplyr::filter(covar_levs == lyr_nm) %>%
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
                                       tmp_terra <- tmp_terra_ls %>%
                                         unlist() %>%
                                         c() %>%
                                         terra::rast() %>%
                                         sum()
                                       
                                       names(tmp_terra) <- lyr_of_interest                    
                                       
                                       return(tmp_terra)
                                     })
      
      data_pred <- tmp_pred_rast_ls %>%
        unlist() %>%
        c() %>%
        terra::rast()
      # %%%%%
      
    }
    
    # otherwise, just combine all into one raster stack
    if (n_prop_covar_lyrs == 0){
      
      tmp_pred_rast <- c(tmp_pred_rast, tmp_bin_rast)[[1]]
      
      # Check all covars in new SpatRaster
      assertthat::assert_that(all(fix_vars %in% names(tmp_pred_rast)))
      # Predict
      if (inherits(obj, "hib_lme4_model")){
        data_pred <- terra_proj_lme4_pred(obj$model, tmp_pred_rast, backtransform)
      }
      if (inherits(obj, "hib_brms_model")){
        data_pred <- terra_proj_brms_pred(obj$model, tmp_pred_rast, backtransform,
                                          incl_obs_uncert)
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
#' @param mod_data A [`data.frame`], of the data used to fit the model
#' @param input_type A [`character`], either "binary" or "non-binary", indicating 
#' the type of data in the raster - i.e., is each grid cell only a single level of 
#' a covariate, or can it contain proportions of multiple levles
#' @keywords data, preparation, model, projection
#' @examples 
#' \dontrun{
#'   r <- terra_covar_prep(covar_nm, covar_df, data_rast, mod_data, input_type = "binary")
#' } 
#' @export
terra_covar_prep <- function(covar_nm, covar_df, data_rast, mod_data, input_type = "binary"){
  # TODO: maybe change to @NoRd
  # TODO: better term than binary
  
  # Assertions
  assertthat::assert_that(
    is.character(covar_nm) || is.factor(covar_nm), 
    is.data.frame(covar_df), 
    inherits(data_rast, "SpatRaster"), 
    is.data.frame(mod_data),
    is.character(input_type)
  )
  
  input_type <- match.arg(input_type, c("binary", "non-binary"), several.ok = FALSE)
  
  # Subset to covar of interest
  df <- covar_df %>%
    dplyr::filter(covar_nm == covar_nm &
                    !is.na(covar_levs_rast))
  
  r <- data_rast[[df$covar_levs_rast]]
  
  if (input_type == "binary"){
    if (max(terra::values(sum(r))) != 1){
      message(crayon::yellow("Warning: Binary SpatRaster contains cells with >1 category! Layers will be discarded ...\n"))
      return (NULL)
    }  
  }
  
  # get variable level IDs
  fac_lev_df <- data.frame(ID = as.numeric(mod_data[[covar_nm]]),
                           category = mod_data[[covar_nm]]
  ) %>%
    unique()
  
  r <- lapply(r, function(lyr, covar_df = df, lev_df = fac_lev_df){
    
    # Need to convert values > 0 to the number that  matches that factor level 
    # in the model df; e.g. primin = 1 in model levels, needs to be 1 here too 
    covar_level <- covar_df %>%
      dplyr::filter(covar_levs_rast == names(lyr)) %>%
      dplyr::pull(covar_levs)
    
    covar_level_numeric <- lev_df$ID[lev_df$category == covar_level]
    
    terra::values(lyr)[terra::values(lyr) > 0] <- covar_level_numeric
    
    # TODO: is ".," needed in terra::categories ??
    if (input_type == "non-binary"){
      lyr <- lyr %>%
        terra::as.factor() %>%
        terra::categories(., value = lev_df)
      # %%%
    }
    return(lyr)
  })
  
  r <- unlist(r) %>%
    terra::rast() 
  
  if (input_type == "binary"){
    r <- r %>%
      max() %>%
      terra::as.factor() %>%
      terra::categories(.,
                        layer = 0,
                        value = fac_lev_df)
  }
  
  # %%%
  names(r) <- rep(covar_nm, dim(r)[3]) 
  return(r)
}


#' # A helper function to project a model using SpatRaster data
#' 
#' @description 
#' This function projects/predicts the supplied model using the SpatRaster 
#' @details 
#' Expects a hib lme4 model and a [`SpatRaster`].
#' @param mod A fitted hib_lme4_mod object 
#' @param data A [`SpatRaster`], prepared for projection
#' @param backtransform A [`function`] that can backtransform the response variable
#' to its original/natual scale.
#' @keywords model, projection
#' @examples 
#' \dontrun{
#'   pred <- terra_proj_lme4_pred(mod, data, backtransform)
#' } 
#' @export
terra_proj_lme4_pred <- function(mod, data, backtransform){
  # TODO: maybe change to @NoRd
  # Assertions
  assertthat::assert_that(
    inherits(mod, "merMod"),
    inherits(data, "SpatRaster"), 
    is.null(backtransform) || is.function(backtransform)
  )
  # Check all fixef covars are in data - already done before func call!
  
  # %%%
  # Project
  # TODO: could also incl df.. - more general
  
  data_pred <- terra::predict(data,
                              mod,
                              re.form = NA)
  names(data_pred) <- "prediction"
  
  # %%%%
  tmp_rast_df <- terra::as.data.frame(data, xy = T)
  
  boot <- lme4::bootMer(mod, 
                        function(m, rast = data){
                          tmp <- terra::predict(rast, m, re.form = NA) %>%
                            terra::as.data.frame()
                          tmp[,1]
                        }, 
                        nsim = 100)
  
  prediction_sd <- apply(boot$t, 2, sd)
  prediction_median <- apply(boot$t, 2, function(x) as.numeric(quantile(x, probs = 0.5, na.rm = TRUE)))
  prediction_q2.5 <- apply(boot$t, 2, function(x) as.numeric(quantile(x, probs = 0.025, na.rm = TRUE)))
  prediction_q97.5 <- apply(boot$t, 2, function(x) as.numeric(quantile(x, probs= 0.975, na.rm = TRUE)))
  
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
  
  # %%%
  if (!is.null(backtransform)){
    data_pred <- backtransform(data_pred)
  }
  # %%% 
  return(data_pred)
}



# %%%%
#' # A helper function to project a model using SpatRaster data
#' 
#' @description 
#' This function projects/predicts the supplied model using the SpatRaster 
#' @details 
#' Expects a hib lme4 model and a [`SpatRaster`].
#' @param mod A fitted hib_brms_mod object 
#' @param data A [`SpatRaster`], prepared for projection
#' @param backtransform A [`function`] that can backtransform the response variable
#' to its original/natual scale.
#' @param incl_obs_uncert A [`logical`] indicating whether or not observation-level
#' uncertainty should be included in the predictions. Setting this parameter to TRUE
#' will produce predictions with larger uncertainty.
#' @keywords model, projection
#' @examples 
#' \dontrun{
#'   pred <- terra_proj_brms_pred(mod, data, backtransform)
#' } 
#' @export
terra_proj_brms_pred <- function(mod, data, backtransform, incl_obs_uncert){
  # TODO: maybe change to @NoRd
  
  # Assertions
  assertthat::assert_that(
    inherits(mod, "brmsfit"),
    inherits(data, "SpatRaster"), 
    is.null(backtransform) || is.function(backtransform),
    is.logical(incl_obs_uncert)
  )
  # Check type of inputs
  # Check all fixef covars are in data 
  
  # %%%
  # Project
  # TODO: could also incl df.. - more general
  
  
  # TODO: 
  # include median for prediction (incl_obs_uncert)
  # %%%%%
  # %%%%%
  if (incl_obs_uncert == T){
    mean_pred <- terra::predict(data,
                                mod,
                                re_form = NA)
    med_pred <- terra::predict(data,
                               mod,
                               re_form = NA,
                               robust = T)
    
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
  if (incl_obs_uncert == F){
    tmp_df <- terra::as.data.frame(data, xy = T)
    mean_pred <- fitted(mod, tmp_df, 
                        re_form = NA) %>%
      as.data.frame()
    median_pred <- fitted(mod, tmp_df, 
                          re_form = NA,
                          robust = T) %>%
      as.data.frame()
    
    
    mean_pred <- mean_pred %>%
      dplyr::rename(prediction_mean  = Estimate,
                    prediction_sd = Est.Error,
                    prediction_q2.5  = Q2.5,
                    prediction_q97.5 = Q97.5)
    median_pred <- median_pred %>%
      dplyr::rename(prediction_median  = Estimate,
                    prediction_mad = Est.Error)
    tmp_pred <- dplyr::bind_cols(tmp_df[,c("x", "y")], 
                                 mean_pred,
                                 median_pred %>%
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
#' @keywords model, projections
#' @examples 
#' \dontrun{
#'   tmp_df <- terra_to_df_inla(obj, data, fix_vars)
#' } 
#' @export
terra_to_df_inla <- function(data, obj, fix_vars){
  # TODO: maybe change to @NoRd
  # TODO: maybe co-opt to use fitted for brms rather then predict - for terra input
  
  # assertions
  assertthat::assert_that(
    inherits(obj, "hib_inla_model"),
    inherits(data, "SpatRaster"), 
    is.character(fix_vars)
  )
  # Covar checks
  
  if (all(fix_vars %in% names(data))){
    # convert terra to df...
    data_df <- terra::as.data.frame(data, xy = T)
  } else {
    # Any layer with a name in fix_vars is good to use 
    covars_in_rast <- names(data)[names(data) %in% fix_vars]
    tmp_pred_rast <- NULL
    n_pred_rast_lyrs <- 0
    if (length(covars_in_rast) > 0 ){
      tmp_pred_rast <- data[[covars_in_rast]]  
      n_pred_rast_lyrs <- dim(tmp_pred_rast)[3]
    }
    
    # any layer with a name not in fix_covars needs to be aligned to a var
    non_covars_in_rast <- names(data)[!names(data) %in% fix_vars]
    tmp_rast <- NULL
    if (length(non_covars_in_rast) > 0 ){
      tmp_rast <- data[[non_covars_in_rast]]
    }
    
    # Assume that layer names match levels of covariates
    # e.g "LUIPrimary_minimal use" or just "Primary_minimal use"
    covars_not_in_rast <- fix_vars[!fix_vars %in% names(data)]
    
    if (length(non_covars_in_rast) > 0 ){
      # for each covar, get unique levels (if factor in mod_df)
      covars_not_in_rast <- purrr::map_dfr(covars_not_in_rast,
                                           function(x){
                                             if (is.factor(obj$data[[x]])){
                                               tmp_out <- data.frame(covar_nm = x,
                                                                     covar_levs = levels(obj$data[[x]]))
                                               # Could add covar_lev_idx here... - save doing later...
                                             }
                                             else {
                                               tmp_out <- data.frame(covar_nm = x,
                                                                     covar_levs = NA)
                                             }
                                             tmp_out
                                           })
      
      covars_not_in_rast <- covars_not_in_rast %>%
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
                                                     is_binary_layer = all(unique(terra::values(tmp_rast[[nm]])) %in% c(0,1))
                                          )
                                        })
      
      covars_not_in_rast <- covars_not_in_rast %>%
        dplyr::left_join(is_binary_layer) %>%
        dplyr::group_by(covar_nm) %>%
        dplyr::mutate(is_binary_covar = all(is_binary_layer, na.rm = T))
      
      # Split binary and non binary layers
      bin_covars <- covars_not_in_rast %>%
        dplyr::filter(is_binary_covar == T)
      
      non_bin_covars <- covars_not_in_rast %>%
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
                                   mod_data = obj$data,
                                   input_type = "binary") %>%
          unlist() %>%
          c() %>%
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
                                    mod_data = obj$data,
                                    input_type = "non-binary") %>%
          unlist() %>%
          c() %>%
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
      
      non_bin_sub <- non_bin_covars %>%
        dplyr::filter(!is.na(covar_levs_rast))
      
      non_bin_r <- tmp_rast[[non_bin_sub$covar_levs_rast]]
      
      tmp_pred_df_ls <- vector("list", n_prop_covar_lyrs)
      names(tmp_pred_df_ls) <- rep("", n_prop_covar_lyrs)
      for (i in 1:n_prop_covar_lyrs){
        # TODO: need to be careful combining these, if first item is null,
        # output becomes a list of 1...
        tmp_r <- c(tmp_pred_rast, tmp_bin_rast, tmp_prop_rast[[i]])[[1]]
        # Convert to df
        tmp_df <- terra::as.data.frame(tmp_r,
                                       xy = T)
        
        # Get proportional weighting
        lyr_nm <- terra::cats(tmp_prop_rast)[[i]] %>% 
          dplyr::filter(ID == as.vector(unname(unique(terra::values(tmp_prop_rast[[i]])))) ) %>% 
          dplyr::pull(category)
        
        lyr_orig_nm <- non_bin_sub %>%
          dplyr::filter(covar_levs == lyr_nm) %>%
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
  
  return(data_df)  
  
}


# TODO:
# Sub functions:
# hib_mod_proj_lme4
# hib_mod_proj_inla



