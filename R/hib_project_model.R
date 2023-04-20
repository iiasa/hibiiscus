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
    is.data.frame(data) || inherits(data, "SpatRaster"), # or a list of SpatRasters
    is.function(backtransform) || is.null(backtransform),
    is.logical(rescale_pred)
  )
  
  
  # Check type of model
  # for now, focus on lmer
  if (inherits(obj, "hib_lme4_model")){
    message(crayon::green("lme4 model detected! Projection in progress...\n"))
    # TODO: wrap in hib_mod_proj_lme4
    
    # Get baseline value for optional rescaling
    # get model intercept
    base_val <- obj$fixedef_coefficients$Coefficient_value[obj$fixedef_coefficients$Coefficient_name == "(Intercept)"]
    # if backtrans - backrans intercept to get baseline val
    if (!is.null(backtransform)){
      base_val <- backtransform(base_val)
    }
    

    # Check type of data 
    if (is.data.frame(data)){
      # Project
      data$prediction <- lme4:::predict.merMod(obj$model,
                                               data,
                                               re.form = NA)
      if (!is.null(backtransform)){
        data$prediction <- backtransform(data$prediction)
      }
      if (rescale_pred){
        # divide predictions by baseline val
        data$prediction_rescaled <- data$prediction / base_val
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
      
      # TODO: what about temporal projections - e.g. over X years...
      # assume data is a list of SpatRaster items, and apply projection process to each...
      # Need to wrap all of raster-focused methods and apply per list item...
      
      # project_terra <- function(obj, 
      #                           data, 
      #                           backtransform, 
      #                           rescale_pred){
      #   
      #   return(data)
      # }
      
      # check if all coef names from model are in raster names -
      # it they are, use terra::predict
      all_covars <- all.vars(formula(delete.response(terms(formula(obj$model)))))
      
      # But what about random effects... - 
      rand_vars <- names(lme4::ranef(obj$model))
      fix_vars <- all_covars[!all_covars %in% rand_vars]
      
      # For initial dev, only interested in fixed effects...
      
      if (all(fix_vars %in% names(data))){
        data_pred <- terra::predict(data,
                                obj$model,
                                re.form = NA)
        
        if (!is.null(backtransform)){
          data_pred <- backtransform(data_pred)
        }
        names(data_pred) <- "prediction"
        
      } else {
        # if not, terra::predict won't work!!
        # Need an alternative approach - custom...
        # message(crayon::yellow("Warning: SpatRaster with non-covariate layers detected! Projections are not yet supported...\n"))
        # stop("In development")
        # TODO: and print statements to help debugging
        
        # Any layer with a name in all_covars is good to use 
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
                             # Could add covar_lev_idx here... - save doing later...
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

          # if multiple covars, but, each layer is binary, that's fine, issue is if
          # multiple covars have proportional cover layers...
          # Check for binary layers...
          is_binary_layer <- purrr::map_dfr(names(tmp_rast),
                                            function(nm){
                                              data.frame(covar_levs_rast = nm,
                                                         is_binary_layer = all(unique(values(tmp_rast[[nm]])) %in% c(0,1))
                                                         )
                                            })
          
          covars_not_in_rast <- covars_not_in_rast %>%
            dplyr::left_join(is_binary_layer) 

          
          covars_not_in_rast <- covars_not_in_rast %>%
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
                                       function(x, 
                                                df = bin_covars,
                                                r = tmp_rast){
                                         df_sub <- df %>%
                                           dplyr::filter(covar_nm == x &
                                                           !is.na(covar_levs_rast))
                                         # print(df_sub)
                                         # print(r)
                                         r_sub <- r[[df_sub$covar_levs_rast]]
                                         # print(r_sub)
                                         if (max(values(sum(r_sub))) != 1){
                                           message(crayon::yellow("Warning: Binary SpatRaster contains cells with >1 category! Layers will be discarded ...\n"))
                                           return (NULL)
                                         }
                                         # convert each layer to factor
                                         fac_lev_df <- data.frame(ID = as.numeric(obj$data[[x]]),
                                                                  category = obj$data[[x]]
                                         ) %>%
                                           unique()
                                         
                                         r_sub <- lapply(r_sub, function(lyr, covar_df = df_sub, lev_df = fac_lev_df){
                                           
                                           # Need to convert binary values (1) to the number that 
                                           # matches that factor level in the model df 
                                           # e.g. primin = 1 in model levels, needs to be 1 here too 
                                           covar_level <- covar_df %>%
                                             dplyr::filter(covar_levs_rast == names(lyr)) %>%
                                             pull(covar_levs)
                                           
                                           covar_level_numeric <- lev_df$ID[lev_df$category == covar_level]
                                           
                                           values(lyr)[values(lyr) > 0] <- covar_level_numeric
                                           
                                           return(lyr)
                                         })
                                         
                                         r_sub <- unlist(r_sub) %>%
                                           terra::rast() %>%
                                           max() %>%
                                           terra::as.factor() %>%
                                           terra::categories(., value = fac_lev_df)
                                         
                                         names(r_sub) <- x
                                         
                                         return(r_sub)
                                       }) %>%
              unlist() %>%
              c() %>%
              terra::rast()
            n_bin_covar_layers <- dim(tmp_bin_rast)[3]
          }
          
          # print(tmp_bin_rast)
  
          
          # For non-binary layers - factor is layer name, minus covar nm if name is covar_nm_lev
          # BUT, also need to record values for later aggregation
          # tmp_prop_rast <- 
          # Actually, extract as above for each layer - but don't max, keep levels in separate layers/stacks
          # One prop stack, one binary stack - predict using latter (binary), then backtrans/aggregate 
          # based on former (prop, weighted av)
          # currently assume only one non_bin_covar...
          
          # print(non_bin_covars)
          
          tmp_prop_rast <- NULL
          n_prop_covar_lyrs <- 0
          if (nrow(non_bin_covars) > 0){
            # No need for purrr
            tmp_prop_rast <- purrr::map(unique(non_bin_covars$covar_nm),
                                        function(x, 
                                                 df = non_bin_covars,
                                                 r = tmp_rast){
                                          df_sub <- df %>%
                                            dplyr::filter(covar_nm == x &
                                                            !is.na(covar_levs_rast))
                                          
                                          r_sub <- r[[df_sub$covar_levs_rast]]
                                          
                                          # convert each layer to factor
                                          fac_lev_df <- data.frame(ID = as.numeric(obj$data[[x]]),
                                                                   category = obj$data[[x]]
                                          ) %>%
                                            unique()
                                          
                                          r_sub <- lapply(r_sub, function(lyr, covar_df = df_sub, lev_df = fac_lev_df){
                                            
                                            # Need to convert prop values to the number that 
                                            # matches that factor level in the model df 
                                            # e.g. primin = 1 in model levels, needs to be 1 here too 
                                            covar_level <- covar_df %>%
                                              dplyr::filter(covar_levs_rast == names(lyr)) %>%
                                              pull(covar_levs)
                                            
                                            covar_level_numeric <- lev_df$ID[lev_df$category == covar_level]
                                            
                                            values(lyr)[values(lyr) > 0] <- covar_level_numeric
                                            
                                            lyr <- lyr %>%
                                              terra::as.factor() %>%
                                              terra::categories(., value = fac_lev_df)
                                            
                                            return(lyr)
                                          })
                                          
                                          r_sub <- unlist(r_sub) %>%
                                            terra::rast() 
                                          
                                          names(r_sub) <- rep(x, dim(r_sub)[3])
                                          
                                          return(r_sub)
                                        }) %>%
              unlist() %>%
              c() %>%
              terra::rast()
            
            # print(cats(tmp_prop_rast))
            # print(tmp_prop_rast)
            # Want to end up with N raster stacks - 1 per prop level
            n_prop_covar_lyrs <- dim(tmp_prop_rast)[3]
          }
          # Check tmp_prop_rast is a terra SpatRaster
          
          
          
          
          # %%%%%%%%%%%%%%%%%
          # TODO : Think this section can be removed...
          # covar_levs_in_rast <- NULL
          # if (sum(covars_not_in_rast$covar_levs_in_rast, na.rm = T) > 0){
          #   covar_levs_rast <- tmp_rast[[covars_not_in_rast$covar_levs[covars_not_in_rast$covar_levs_in_rast]]]
          #   
          #   # Need to harmonise these layers - extract each layer individually
          #   
          #   # convert to character/factor based on prop > 0
          #   
          #   # Need to record proportions...
          #   
          # 
          # }
          # 
          # covar_nm_levs_in_rast <- NULL
          # if (sum(covars_not_in_rast$covar_nm_levs_in_rast, na.rm = T) > 0){
          #   covar_nm_levs_rast <- tmp_rast[[covars_not_in_rast$covar_nm_levs[covars_not_in_rast$covar_nm_levs_in_rast]]]
          # }
          # %%%%%%%%%%%%%%%%%
         
          
        }
        
        # tmp_pred_rast + tmp_bin_rast + prop_pred_rast
        # check N layers of each
        # if n_prop_covar_lyrs > 1  - list of rasters - each named based on 
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

            lyr_nm <- terra::cats(tmp_prop_rast)[[i]] %>% 
              dplyr::filter(ID == as.vector(unname(unique(values(tmp_prop_rast[[i]])))) ) %>% 
              pull(category)
            
            names(tmp_pred_rast_ls)[i] <- lyr_nm
            
            tmp_pred_rast_ls[[i]] <- terra::predict(tmp_pred_rast_ls[[i]],
                                                    obj$model,
                                                    re.form = NA)
            if (!is.null(backtransform)){
              tmp_pred_rast_ls[[i]] <- backtransform(tmp_pred_rast_ls[[i]])
            }
            
            lyr_orig_nm <- non_bin_sub %>%
              dplyr::filter(covar_levs == lyr_nm) %>%
              pull(covar_levs_rast)
            
            tmp_pred_rast_ls[[i]] <- non_bin_r[[lyr_orig_nm]] * tmp_pred_rast_ls[[i]]
          }
          
          
          # then aggr - 
          # link original non-binary raster layers to outputs 
          # Output list names match covar levels, 
          # for each list item, get name, get non_bin_covars$covar_levs_rast[non_bin_covars$covar_levs == nm]
          # Original raster layer names may not - these match non_bin_covars$covar_levs_rast
          
          # for each item in list, multiply prediction by original raster (correctly matched names)
          # then sum all raster in list to get weighted av pred...
          data_pred <- tmp_pred_rast_ls %>%
            unlist() %>%
            c() %>%
            terra::rast() %>%
            sum()
          
          names(data_pred) <- "prediction"
        }
        
        # otherwise, just combine all into one
        if (n_prop_covar_lyrs == 0){

          tmp_pred_rast <- c(tmp_pred_rast, tmp_bin_rast)[[1]]

          # Predict
          data_pred <- terra::predict(tmp_pred_rast,
                                 obj$model,
                                 re.form = NA)
          
          if (!is.null(backtransform)){
            data_pred <- backtransform(data_pred)
          }
          names(data_pred) <- "prediction"
        }
        
        # TODO:
        # rescale predictions compared to baseline...
        
        # for any of these - retain layers - drop others?
        
        
        
        # need to extract covar numeric values (LU proportion)
        
      }
      
      data_pred_rescaled <- NULL
      if (rescale_pred){
        # divide predictions by baseline val
        data_pred_rescaled <- data_pred / base_val
        names(data_pred_rescaled) <- "prediction_rescaled"
      }
      
      data <- c(data, data_pred, data_pred_rescaled)
      
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
  # return(data)
  
  # TODO:
  # Return a projection object
  # create projection object - need model obj, projection output, and projection_specs (backrans, rescale)
  message("Creating projection object...\n")
  proj <- hib_mod_projection$new(obj = obj,
                                 proj = data,
                                 backtransform = backtransform,
                                 rescale_pred = rescale_pred)
  return(proj)
}


# TODO:
# Sub functions:
# hib_mod_proj_lme4
# hib_mod_proj_inla

# TODO:
# uncertainty in predictions - easy for INLA - sampling? for lme4
