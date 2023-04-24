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
#' @keywords model, projections
#' @examples 
#' \dontrun{
#'   proj <- hib_mod_project(obj, data)
#' } 
#' @export
hib_mod_project <- function(obj, 
                             data,
                             backtransform = NULL,
                             rescale_pred = F){
  # Assertions 
  # is data a df/tibble/ raster/stack 
  assertthat::assert_that(
    inherits(obj, "hib_lme4_model") || inherits(obj, "hib_inla_model"), # obj is a hib model
    is.data.frame(data) || inherits(data, "SpatRaster") || (is.list(data) & !is.data.frame(data)),
    is.function(backtransform) || is.null(backtransform),
    is.logical(rescale_pred)
  )
  
  
  # Check type of model
  # for now, focus on lmer
  if (inherits(obj, "hib_lme4_model")){
    message(crayon::green("lme4 model detected! Projection in progress...\n"))
    # TODO: wrap in hib_mod_proj_lme4??
    
    # Get baseline value for optional rescaling
    # get model intercept
    base_val <- obj$fixedef_coefficients$Coefficient_value[obj$fixedef_coefficients$Coefficient_name == "(Intercept)"]
    # if backtrans - backrans intercept to get baseline val
    if (!is.null(backtransform)){
      base_val <- backtransform(base_val)
    }
    
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
      
      data$prediction <- lme4:::predict.merMod(obj$model,
                                               data,
                                               re.form = NA)
      # What if we want to include  random effect??
      
      if (!is.null(backtransform)){
        data$prediction <- backtransform(data$prediction)
      }
      if (rescale_pred){
        data$prediction_rescaled <- data$prediction / base_val
      }
    }
    
    # Single raster object
    if (inherits(data, "SpatRaster")){
      data <- terra_proj(data, obj, backtransform, rescale_pred, base_val, fix_vars)
    }
    # List of raster objects
    if ((is.list(data) & !is.data.frame(data))){
      # Apply individual terra proj to each item of list...
      data <- purrr::map(data,
                         terra_proj,
                         obj = obj, 
                                       backtransform = backtransform, 
                                       rescale_pred = rescale_pred, 
                                       base_val = base_val, 
                                       fix_vars = fix_vars)
      # print(data)
      # print(is.list(data))
      # print(!is.data.frame(data))
    }
    
  } # TODO: add inla support
  else if (inherits(obj, "hib_inla_model")){
    message(crayon::yellow("Warning: INLA model detected! Projections are not yet supported...\n"))
    stop("In development")
    
    # TODO: put into hib_mod_proj_inla
 
    
    # Check type of data 
    
    # For df, just incorporate into model re-fit... 
    # - need to know original model formula...
    # need to check presence of fixef, and lat/long if spatial ranef...
    
    # If terra, convert rasters to dfs
    # If list of terra, convert each raster to df, add name of item in list, bind all
    
    # get projection A using mesh...... - and mesh specs... - e.g. 2d v 3d
    
    # new, pred INLA-stack
    # re-fit
    # Project
    
    # if input was terra/list of terra, convert pred df back to terra format
  }
  
  # create projection object - need model obj, projection output, and projection_specs (backrans, rescale)
  message("Creating projection object...\n")
  proj <- hib_mod_projection$new(obj = obj,
                                 proj = data,
                                 backtransform = backtransform,
                                 rescale_pred = rescale_pred)
  return(proj)
}


#' # A helper function to wrap the preparation of SpatRaster data, and then project 
#' a hib, mixed-effects model - currently only for lme4 style model..
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
#' @keywords model, projections
#' @examples 
#' \dontrun{
#'   tmp_proj <- terra_proj(obj, data, backtransform, rescale_pred, base_val, fix_vars)
#' } 
#' @export
terra_proj <- function(data, obj, backtransform, rescale_pred, base_val, fix_vars){
  # assertions
  assertthat::assert_that(
    inherits(obj, "hib_lme4_model") || inherits(obj, "hib_inla_model"),
    inherits(data, "SpatRaster"), 
    is.function(backtransform) || is.null(backtransform),
    is.logical(rescale_pred)
  )
  # Covar checks
  
  if (all(fix_vars %in% names(data))){
    data_pred <- terra::predict(data,
                                obj$model,
                                re.form = NA)
    
    if (!is.null(backtransform)){
      data_pred <- backtransform(data_pred)
    }
    names(data_pred) <- "prediction"
  } else {
    # Any layer with a name in fix_vars is good to use 
    covars_in_rast <- names(data)[names(data) %in% fix_vars]
    tmp_pred_rast <- NULL
    n_pred_rast_lyrs <- 0
    if (length(covars_in_rast) > 0 ){
      tmp_pred_rast <- data[[covars_in_rast]]  
      n_pred_rast_lyrs <- dim(tmp_pred_rast)[3]
    }
    
    # any layer with a name not in all_covars needs to be aligned to a var
    non_covars_in_rast <- names(data)[!names(data) %in% fix_vars]
    tmp_rast <- NULL
    if (length(non_covars_in_rast) > 0 ){
      tmp_rast <- data[[non_covars_in_rast]]
    }
    
    # Assume that layer names match levels of covariates
    # e.g "LUIPrimary_minimal" use or just "Primary_minimal use"
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
                      covar_levs_rast = case_when(covar_levs_in_rast == T ~ covar_levs,
                                                  covar_nm_levs_in_rast == T ~ covar_nm_levs,
                                                  T ~ NA))
      
      # Check for binary layers...
      is_binary_layer <- purrr::map_dfr(names(tmp_rast),
                                        function(nm){
                                          data.frame(covar_levs_rast = nm,
                                                     is_binary_layer = all(unique(values(tmp_rast[[nm]])) %in% c(0,1))
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
        
        tmp_pred_rast_ls[[i]] <- terra_proj_lme4_pred(obj$model, 
                                                      tmp_pred_rast_ls[[i]], 
                                                      backtransform)
        
        # Get proportional weighting
        lyr_nm <- terra::cats(tmp_prop_rast)[[i]] %>% 
          dplyr::filter(ID == as.vector(unname(unique(values(tmp_prop_rast[[i]])))) ) %>% 
          pull(category)
        
        lyr_orig_nm <- non_bin_sub %>%
          dplyr::filter(covar_levs == lyr_nm) %>%
          pull(covar_levs_rast)
        
        # multiply prediction by original raster (correctly matched names)
        tmp_pred_rast_ls[[i]] <- non_bin_r[[lyr_orig_nm]] * tmp_pred_rast_ls[[i]]
        
        names(tmp_pred_rast_ls)[i] <- lyr_nm
        
      }
      
      # then sum all raster in list to get weighted av pred...
      data_pred <- tmp_pred_rast_ls %>%
        unlist() %>%
        c() %>%
        terra::rast() %>%
        sum()
      
      names(data_pred) <- "prediction"
    }
    
    # otherwise, just combine all into one raster stack
    if (n_prop_covar_lyrs == 0){
      
      tmp_pred_rast <- c(tmp_pred_rast, tmp_bin_rast)[[1]]
      
      # Check all covars in new SpatRaster
      assertthat::assert_that(all(fix_vars %in% names(tmp_pred_rast)))
      # Predict
      data_pred <- terra_proj_lme4_pred(obj$model, 
                                        tmp_pred_rast, 
                                        backtransform)
    }
  }
  
  # Optional rescaling
  data_pred_rescaled <- NULL
  if (rescale_pred){
    # divide predictions by baseline val
    data_pred_rescaled <- data_pred / base_val
    names(data_pred_rescaled) <- "prediction_rescaled"
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
  # TODO: better term than binary
  
  # Assertions
  assertthat::assert_that(
    is.character(covar_nm) || is.factor(covar_nm), 
    is.data.frame(covar_df), 
    inherits(data_rast, "SpatRaster"), 
    is.data.frame(mod_data),
    is.character(input_type)
  )
  
  
  # Subset to covar of interest
  df <- covar_df %>%
    dplyr::filter(covar_nm == covar_nm &
                    !is.na(covar_levs_rast))
  
  r <- data_rast[[df$covar_levs_rast]]
  
  if (input_type == "binary"){
    if (max(values(sum(r))) != 1){
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
      pull(covar_levs)
    
    covar_level_numeric <- lev_df$ID[lev_df$category == covar_level]
    
    values(lyr)[values(lyr) > 0] <- covar_level_numeric
    
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
#' to it original/natual scale.
#' @keywords model, projection
#' @examples 
#' \dontrun{
#'   pred <- terra_proj_lme4_pred(mod, data, backtransform)
#' } 
#' @export
terra_proj_lme4_pred <- function(mod, data, backtransform){
  # Assertions
  # Check type of inputs
  # Check all fixef covars are in data 
  
  # %%%
  # Project
  # TODO: could also incl df.. - more general
  
  data_pred <- terra::predict(data,
                              mod,
                              re.form = NA)
  if (!is.null(backtransform)){
    data_pred <- backtransform(data_pred)
  }
  names(data_pred) <- "prediction"
  # %%% 
  return(data_pred)
}
 


# terra_proj_inla_pred <- function(obj, data, backtransform, rescal_pred){
#   # Assertions...
#   
# }   

# TODO:
# uncertainty in predictions - easy for INLA - sampling? for lme4

# TODO:
# Sub functions:
# hib_mod_proj_lme4
# hib_mod_proj_inla



