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
#' @keywords model, projections
#' @examples 
#' \dontrun{
#'   proj <- hib_mod_project(obj, data)
#' } 
#' @export
hib_mod_project <- function(obj, 
                            data,
                            # TODO: take in function
                            backtransform = NULL){
  # Assertions 
  # is data a df/tibble/ raster/stack 
  assertthat::assert_that(
    inherits(obj, "hib_lme4_model") || inherits(obj, "hib_inla_model"), # obj is a hib model
    is.data.frame(data) || inherits(data, "SpatRaster"),
    is.function(backtransform) || is.null(backtransform)
  )
  
  
  # Check type of model
  # for now, focus on lmer
  if (inherits(obj, "hib_lme4_model")){
    message(crayon::green("lme4 model detected! Projection in progress...\n"))
    # TODO: wrap in hib_mod_proj_lme4
    
    # Check type of data 
    if (is.data.frame(data)){
      # Project
      data$prediction <- lme4:::predict.merMod(obj$model,
                                               data,
                                               re.form = NA)
      if (!is.null(backtransform)){
        data$prediction <- backtransform(data$prediction)
      }
    }
    # TODO:
    # Above assumes each row is a unique "site" with single LUI + other vars..
    # What if each site is a grid cell with multiple LUIs as different proportions
    # 
    # Current approach requires df input...
    
    
    if (inherits(data, "SpatRaster")){
      # message(crayon::yellow("Warning: SpatRaster input detected! Projection are not yet supported...\n"))
      # stop("In development")
      
      # TODO: check if all coef names from model are in raster names -
      # it they are, use terra::predict
      all_covars <- all.vars(formula(delete.response(terms(formula(obj$model)))))
      
      # But what about random effects... - 
      rand_vars <- names(lme4::ranef(obj$model))
      fix_vars <- all_covars[!all_covars %in% rand_vars]
      
      # For initial dev, only interested in fixed effects...
      if (all(fix_vars %in% names(data))){
        data <- terra::predict(data,
                                obj$model,
                                re.form = NA)
        
        if (!is.null(backtransform)){
          data <- backtransform(data)
        }
      } else {
        # if not, terra::predict won't work!!
        # Need an alternative approach - custom...
        message(crayon::yellow("Warning: SpatRaster with non-covariate layers detected! Projections are not yet supported...\n"))
        stop("In development")
        
        # Any layer with a name in all_covars is good to use 
        covars_in_rast <- names(data)[names(data) %in% fix_vars]
        tmp_pred_rast <- NULL
        if (length(covars_in_rast) > 0 ){
          tmp_pred_rast <- data[[covars_in_rast]]  
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
        
        # if (length(covars_not_in_rast) > 1){
        #   message(crayon::yellow("Warning: SpatRaster with non-covariate layers detected! Projections are not yet supported...\n"))
        #   stop("In development")
        #   
        # }
        # What if covars_not_in_rast is > 1 - i.e. multiple covariates are represented
        # by multiple raster layers (assuming proportional cover for each covar level in
        # a distinct raster layer) - how to predict if specific combinations of
        # covar levels and proportions are not known??
        # Check it mathematically...simple 0.5, 0.5 per covar, but one is vertical
        # split and other is horizontal, v both vertical/horizontal - 
        # Have to assume that this wont be the case??
        
        
        # what if multiple factor covars - with different proportions... how to aggregate?!
        # e.g., LU is 0.2 Primary, 0.8, Cropland, and UI is 0.4 Minimal, 0.6 Light
        # Need to know proportional cover of combined 
        # e.g. 0.2 Primary Minimal, 0.2 Cropland Minimal, 0.6 Cropland Light) 
        # - otherwise impacts can't be properly estimated??.
        
        if (length(non_covars_in_rast) > 0 ){
          # for each covar, get unique levels (if factor in mod_df)
          covars_not_in_rast <- purrr::map_dfr(covars_not_in_rast,
                         function(x){
                           if (is.factor(obj$data[[x]])){
                             tmp_out <- data.frame(covar_nm = x,
                                                   covar_levs = levels(obj$data[[x]]))
                           }
                           else {
                             tmp_out <- data.frame(covar_nm = x,
                                                   covar_levs = NA)
                           }
                           tmp_out
                         })
          # Combine covar nm with levels
          covars_not_in_rast$covar_nm_levs <- paste0(covars_not_in_rast$covar_nm,
                                                     covars_not_in_rast$covar_levs)
          # check for presence in tmp_rast
          covars_not_in_rast$covar_levs_in_rast <- covars_not_in_rast$covar_levs %in% names(tmp_rast)
          covars_not_in_rast$covar_nm_levs_in_rast <- covars_not_in_rast$covar_nm_levs %in% names(tmp_rast)
          
          covars_not_in_rast <- covars_not_in_rast %>%
            dplyr::mutate(covar_levs_rast = case_when(covar_levs_in_rast == T ~ covar_levs,
                                                      covar_nm_levs_in_rast == T ~ covar_nm_levs,
                                                      T ~ NA))
          un_covars_not_in_rast <- covars_not_in_rast %>% 
            dplyr::filter(covar_levs_in_rast == T |
                            covar_nm_levs_in_rast) %>%
            dplyr::pull(covar_nm) %>%
            unique()
          
          # if multiple covars, but, each layer is binary, that's fine, issue is if
          # multiple covars have proportional cover layers...
          # Check for binary layers...
          is_binary_layer <- purrr::map(names(tmp_rast),
                                        function(nm){
                                          all(unique(values(tmp_rast[[nm]])) %in% c(0,1))
                                          # all( %in% c(0, 1))
                                        }) %>%
            unlist()
            
          is_binary_layer <- data.frame(covar_levs_rast = names(tmp_rast),
                                        is_binary_layer = is_binary_layer)
          
          covars_not_in_rast <- covars_not_in_rast %>%
            dplyr::left_join(is_binary_layer) 
          
          bin_non_bin_covars <- covars_not_in_rast %>%
            dplyr::group_by(covar_nm) %>%
            dplyr::summarise(is_binary_covar = all(is_binary_layer, na.rm = T))
          
          covars_not_in_rast <- covars_not_in_rast %>%
            dplyr::left_join(bin_non_bin_covars)
          
          # Split binary and non binary layers
          bin_covars <- covars_not_in_rast %>%
            dplyr::filter(is_binary_covar == T)
          
          non_bin_covars <- bin_non_bin_covars %>%
            dplyr::filter(is_binary_covar == F)
          
          if (length(non_bin_covars) > 1){
            message(crayon::yellow("Warning: SpatRaster with multiple non-binary layers detected! Projections are not yet supported...\n"))
            stop("In development")
          }
          
          
          # For binary layers - factor is layer name, minus covar nm if name is covar_nm_lev
          # Collapse to single stack with layer names as covar_nm
          modal
          
          purrr::map(unique(bin_covars$covar_nm),
                     function(x, 
                              df = bin_covars,
                              r = tmp_rast){
                       df_sub <- df %>%
                         dplyr::filter(covar_nm == x &
                                         !is.na(covar_levs_rast))
                       # print(head(df_sub))
                       r_sub <- r[[df_sub$covar_levs_rast]]
                       if (max(values(sum(r_sub))) != 1){
                         message(crayon::yellow("Warning: Binary SpatRaster contains cells with >1 category! Layers will be discarded ...\n"))
                         return (NULL)
                       }
                       # convert each layer to factor
                       r_sub <- lapply(r_sub, function(lyr){
                         # print(names(lyr))
                         tmp_r <- terra::as.factor(lyr)
                         tmp_r <- categories(tmp_r, 
                                    value = data.frame(ID = 1,
                                                       category = names(lyr)))
                         names(tmp_r) <- names(lyr)
                         
                         return(tmp_r)
                       })
                       return(r_sub)
                     })
          
          # For non-binary layers - factor is layer name, minus covar nm if name is covar_nm_lev
          # BUT, also need to record values for later aggregation
          
          # Want to end up with N raster stacks - 1 per prop level
          # tmp_pred_rast + binary_pred_rast + prop_pred_rast
          
          # Then backtrans - if wanted
          
          # then aggr - 
          
          covar_levs_in_rast <- NULL
          if (sum(covars_not_in_rast$covar_levs_in_rast, na.rm = T) > 0){
            covar_levs_rast <- tmp_rast[[covars_not_in_rast$covar_levs[covars_not_in_rast$covar_levs_in_rast]]]
            
            # Need to harmonise these layers - extract each layer individually
            
            # convert to character/factor based on prop > 0
            
            # Need to record proportions...
            

          }
          
          covar_nm_levs_in_rast <- NULL
          if (sum(covars_not_in_rast$covar_nm_levs_in_rast, na.rm = T) > 0){
            covar_nm_levs_rast <- tmp_rast[[covars_not_in_rast$covar_nm_levs[covars_not_in_rast$covar_nm_levs_in_rast]]]
          }
          
          
        }
        
        
        
        
        
        # for any of these - retain layers - drop others?
        
        
        
        # need to extract covar numeric values (LU proportion)
        
      }
      
      
      # If a raster, need to have layers as cat vars (either binary or prop of cell) + cont vars
      # cat layers need to be factors
      # layers corresponding to covars that are untransformed?? - extract as is...
      
      # any layer with a name in all
      # need to convert raster to df?, add prop col and binary col for cat vars
      # make long - i.e., row per categorical var combination... 
      # predict, 
      # if back trans...
      # aggr per cell
      # return raster of input dims
      
    }
    
  } # TODO: add inla support
  else if (inherits(obj, "hib_inla_model")){
    message(crayon::yellow("Warning: INLA model detected! Projections are not yet supported...\n"))
    stop("In development")
    
    # TODO: put into hib_mod_proj_inla
    # Check type of data 
    # Project
  }

  # Return projected data
  return(data)
  
  # TODO:
  # Return a projection object
}


# TODO:
# Sub functions:
# hib_mod_proj_lme4
# hib_mod_proj_inla

# TODO:
# uncertainty in predictions - easy for INLA - sampling? for lme4
