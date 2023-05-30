#' Cross-validate a fitted hib mixed-effects model
#' 
#' @description 
#' This function re-fits models on different subsets of the data, so that subsequent 
#' evaluation of model performance can take place.
#' @details 
#' Expects a hib model object...
#' @param mod_obj A [`hib_x_model`] object.
#' @param db_obj A [`hib_database`] object, or, if modelling used multiple data-
#' sources, a list of database objects
#' @param fold_type A [`character`] indicating the type of cross-validation; one
#' of: ..., <spatial_options>, ..., "Custom" - ?"Random-5", "Random-10"?
#' @param custom_col A [`character`] indicating the column in the model data to 
#' use to define cross-validation folds, when fold_type == "custom"
#' @param backtransform A function to convert the model prediction back to the original
#' scale of interest
#' @keywords model, evaluation
#' @examples 
#' \dontrun{
#'   cv_out <- hib_mod_cv(mod_obj, db_obj, fold_type = "Source_ID")
#' } 
#' @export
hib_mod_cv <- function(mod_obj, db_obj, fold_type, custom_col = NULL,
                       backtransform = NULL){
  # Assertion checks
  assertthat::assert_that(
    inherits(mod_obj, "hib_lme4_model") || inherits(mod_obj, "hib_inla_model") || inherits(mod_obj, "hib_brms_model"),     
    inherits(db_obj, "Database") || ( is.list(db_obj) & all(unlist(lapply(obj, function(x){inherits(x, "Database")}))) ),    
    is.character(fold_type),
    is.character(custom_col) || is.null(custom_col),
    is.function(backtransform) || is.null(backtransform)
  )
  
  # allow multiple fold types - 
  # also allow multiple custom cols - 
  # make a spec df... - each of fold types (unique), plus as many custom as unique
  # custom cols
  
  fold_type <- match.arg(fold_type,
                         c("Random-5", "Random-10",
                           "Source_ID", "SS", #"SSBS",
                           "Country", "UN_subregion", "UN_region",
                           "Biome", "Realm",
                           "Custom"), # loo??
                         several.ok = T)
  fold_type <- unique(fold_type)
  
  if (!is.null(custom_col)){
    if (!"Custom" %in% fold_type){
      message(crayon::red(paste0("Warning: custom columns provided, but 'Custom' ",
                                 "cross-validation not specified. ",
                                 "Custom cross-validation will be run.")))
    }
    fold_type <- c(fold_type, "Custom")
  }
  if ("Custom" %in% fold_type & is.null(custom_col)){
    message(crayon::red(paste0("Warning: 'Custom' cross-validation specified, but",
                               "custom columns not provided. ",
                               "'Custom' cross-validation will NOT be run.")))
    fold_type <- fold_type[fold_type != "Custom"]
  }
  
  
  cv_spec <- data.frame(fold_type = fold_type[fold_type != "Custom"],
                        fold_col = fold_type[fold_type != "Custom"]) %>%
    dplyr::mutate(fold_col = dplyr::case_when(fold_type == "Random-5" ~ "fold_id",
                                              fold_type == "Random-10" ~ "fold_id",
                                              T ~ fold_col))
  
  if ("Custom" %in% fold_type){
    cv_spec <- bind_rows(cv_spec,
                         expand.grid(fold_type = c("Custom"),
                                     fold_col = custom_col,
                                     stringsAsFactors = F) %>%
                           as.data.frame()
                         )
  }
  

  # then with spec sorted, purr map for each spec, and apply relevant func...
  # names as <fold_type>_<custom_col>
  # each list item has pred df, coef df, fold_type, n_folds, and an in-sample pred_df...
  # plus any model-specific specs...
  if (inherits(mod_obj, "hib_lme4_model")){
    
    # TODO: print - cv supplied lme4 model, (formula...)
    message(crayon::green(paste0("Cross-validating supplied lme4 model:\n\t", 
                                 format(stats::formula(mod_obj$model)), "\n")))
    cv_out <- purrr::pmap(cv_spec,
                          function(fold_type, fold_col,  
                                   mod_ = mod_obj,
                                   db_ = db_obj,
                                   backtr = backtransform){
                            hib_cv_lme4(mod_, db_, fold_type, fold_col, backtr)
                          })
    
    # %%%%%
    # In-sample preds
    orig_form <- stats::formula(mod_obj$model)
    insample_df <- mod_obj$data
    insample_df$prediction <- predict(mod_obj$model,
                                      insample_df,
                                      re.form = NA)
    
    # backtransform prediction, y as resp term, 
    if (!is.null(backtransform)){
      insample_df <- insample_df %>%
        dplyr::mutate(prediction = backtransform(prediction),
                      y = !!sym(all.vars(orig_form[[2]])))
    }
    # or convert resp value to y - based on eqn, and don't backtransform prediction
    if (is.null(backtransform)){
      insample_df <- insample_df %>%
        dplyr::mutate(y = !!orig_form[[2]] )
    }
    
    insample_df <- insample_df %>%
      dplyr::mutate(error = y - prediction)
    
    insample_metrics <- calc_cv_metrics(insample_df$prediction, 
                                        insample_df$y, 
                                        insample_df$error)
    
    insample_res = list(insample_pred_df = insample_df, 
                        insample_metrics = insample_metrics)

    # %%%%%
    # cv_out <- hib_cv_lme4(mod_obj, db_obj, fold_type, custom_col,
    #                       backtransform)
  }
  
  if (inherits(mod_obj, "hib_inla_model")){
    stop("INLA cv in development!")
    cv_out <- hib_cv_inla(mod_obj, db_obj, fold_type, custom_col,
                          backtransform)

  }
  
  if (inherits(mod_obj, "hib_brms_model")){
    stop("INLA cv in development!")
    cv_out <- hib_cv_brms(mod_obj, db_obj, fold_type, custom_col,
                          backtransform)
    

  }
  
  # Return a cv object.. - incl original model? + df of cv predictions 
  # initialise cv_output...
  message(crayon::green("\nCollating cross-validation output ..."))
  
  cv_obj <- hib_cv_res$new(
    mod_obj = mod_obj,
    cv_res = cv_out,
    insamp_res = insample_res)
  
  return(cv_obj)
}



#' Cross-validate a fitted hib_lme4_model
#' 
#' @description 
#' This function re-fits the model on different subsets of the data, so that subsequent 
#' evaluation of model performance can take place.
#' @details 
#' Expects a hib_lme4_model object...
#' @param mod_obj A [`hib_lme4_model`] object.
#' @param db_obj A [`hib_database`] object, or, if modelling used multiple data-
#' sources, a list of database objects
#' @param fold_type A [`character`] indicating the type of cross-validation; one
#' of: ..., <spatial_options>, ..., "Custom" - ?"Random-5", "Random-10"?
#' @param custom_col A [`character`] indicating the column in the model data to 
#' use to define cross-validation folds, when fold_type == "custom"
#' @param backtransform A function to convert the model prediction back to the original
#' scale of interest
#' @keywords model, evaluation
#' @examples 
#' \dontrun{
#'   cv_out <- hib_cv_lme4(mod_obj, db_obj, fold_type = "Source_ID")
#' } 
#' @noRd
hib_cv_lme4 <- function(mod_obj, db_obj, fold_type, fold_col,
                        backtransform = NULL){
  
  # Assertion checks 
  assertthat::assert_that(
    inherits(mod_obj, "hib_lme4_model"),
    inherits(db_obj, "Database") || is.list(db_obj), # TODO list of db...
    is.character(fold_type),
    is.character(fold_col),
    is.function(backtransform) || is.null(backtransform)
  )
  
  
  # Check fold_type matches one of the options
  fold_type <- match.arg(fold_type, 
                         c("Random-5", "Random-10", 
                           "Source_ID", "SS", #"SSBS",
                           "Country", "UN_subregion", "UN_region",
                           "Biome", "Realm",
                           "Custom"), # loo??
                         several.ok = F)
  
  # %%%%%
  # TODO: reconfigure to use fold_type and fold_col as supplied
  fold_type_msg <- fold_type
  if (fold_type == "Custom"){
    fold_type_msg <- paste0(fold_type, " (", fold_col, ")")
  }
  message(crayon::green(paste0("\nRunning cross-validation using: ", fold_type_msg, "\n")))
  # fold_col <- NA
  # if (fold_type %in% c("Source_ID", "SS", "Country", "UN_subregion", "UN_region",
  #                      "Biome", "Realm")){
  #   fold_col <- fold_type
  # }
  # if (fold_type == "custom"){
  #   fold_col <- custom_col
  # }
  
  
  # Ensure fold_col is in database, and extract ...
  if (fold_type %in% c("Source_ID", "SS", "Country", "UN_subregion", "UN_region",
                       "Biome", "Realm", "Custom")){
    # Check fold_col is in db_obj$site
    if(inherits(db_obj, "Database")){
      assertthat::assert_that(
        assertthat::has_name(db_obj$sites, fold_col)
      )
      fold_df <- db_obj$sites[,c("SSBS", fold_col)]
    }
    # or each db_obj$site if db_obj is a list
    if(is.list(db_obj)){
      assertthat::assert_that(
        all(unlist(lapply(db_obj, function(x){fold_col %in% colnames(x$sites)})))
      )
      
      fold_df <- lapply(db_obj, function(x){x$sites[,c("SSBS", fold_col)]}) %>%
        dplyr::bind_rows()
    }
  }
  
  # Random cv doesn't seem sensible... - ignores important struture in data - e.g. SS
  # Does that matter??
  if (fold_type == "Random-5"){
    # Randomly sample/shuffle SSBS, assign to 5 indices
    fold_df <- data.frame(SSBS = sample(mod_obj$data$SSBS),
                          fold_id = seq(1,5))
    # fold_col <- "fold_id"
  }

  if (fold_type == "Random-10"){
    # Randomly sample/shuffle SSBS, assign to 10 indices
    fold_df <- data.frame(SSBS = sample(mod_obj$data$SSBS),
                          fold_id = seq(1,10))
    # fold_col <- "fold_id"
  }
  # %%%%%
  
  
  # Should we also try to "keep" SS together in other blockings...
  
  # Need to obtain model fitting dataframe - mod_obj$data
  mod <- mod_obj$model
  # model formula - formula(mod_obj$model)
  orig_form <- stats::formula(mod)
  
  # Other key info from model fitting
  # family # lmer v glmer
  fam <- NULL
  if (is(mod@resp, "glmResp")){
    fam <- mod@resp$family
  }
  
  orig_data <- mod_obj$data
  
  # weights
  incl_weights <- F
  if ("data_weight" %in% colnames(orig_data)){
    incl_weights <- T
  }
  
  # reml v ml
  if (mod@devcomp$dims[["REML"]]){
    REML <- T
  } else {REML <- F}
  
  
  ci_type <- mod_obj$ci_type
  
  # Also need hib database - for wider site data!! 
  # If model is based on contrasts (i.e., comp-sim style), need to get, info on both
  # SSBS1 and SSBS2... - below, am just using SSBS1...
  
  # Add important fold variable - from db_obj - based on SSBS, or, 
  # if SSBS1/SSBS2 in orig_data, use these...
  if ("SSBS" %in% colnames(orig_data)){
    cv_data <- dplyr::left_join(orig_data,
                                fold_df,
                                by = "SSBS")
  }
  if ("SSBS1" %in% colnames(orig_data)){
    cv_data <- dplyr::left_join(orig_data,
                                fold_df,
                                by = c("SSBS1" = "SSBS"))
  }
  
  
  fold_ids <- unique(cv_data %>% dplyr::pull(!!fold_col))
  # print(fold_ids)
  n_folds <- length(fold_ids)
  # Make a list to store output - length = length(unique(new_df[,fold_var]))
  cv_out_ls <- vector("list", n_folds)
  
  coef_out_ls <- vector("list", n_folds)
    
  div_ <- 1
  if (n_folds > 10){
    div_ <- 10
  }
  if (n_folds > 100){
    div_ <- 100
  }
  # For every fold - subset data
  
  message(crayon::green(paste0("\tTotal number of folds: ", n_folds, "\n\n\tProgress")))
  
  pb <- txtProgressBar(min = 1, max = n_folds, style = 3)
  
  for (i in 1:n_folds){
    
    # if (i %% div_ == 0){
    #   message(crayon::green(paste0("Fold ", i, " of ", n_folds, "...\n")))
    # }
    # print(head(cv_data))
    # print(i)
    # print(fold_ids[i])

    # tr_idx <- cv_data %>% dplyr::pull(!!fold_col) != fold_ids[i]
    # te_idx <- cv_data %>% dplyr::pull(!!fold_col) == fold_ids[i]
    # 
    # tr <- cv_data[tr_idx,]
    # te <- cv_data[te_idx,]
    
    tr <- cv_data %>%
      dplyr::filter(!!sym(fold_col) != fold_ids[i])
    te <- cv_data %>%
      dplyr::filter(!!sym(fold_col) == fold_ids[i])
    
    # print(head(tr))
    # print(head(te))
    
    # Make a fold_id/fold_idx column for te...
    te <- te %>%
      dplyr::mutate(fold_id = !!sym(fold_col))
    
    # TODO:
    # Need to adjust weights in tr so sum per data type - dataset_modelling_id
    # sum to same value as in orig data...
    if (incl_weights == T){
      weight_df <- orig_data %>%
        dplyr::group_by(dataset_modelling_id) %>%
        dplyr::summarise(prop_weight = sum(data_weight))
      
      tr <- tr %>%
        select(-data_weight)
      te <- te %>%
        select(-data_weight)
      
      tr <- dplyr::left_join(tr, weight_df) %>%
        dplyr::group_by(dataset_modelling_id) %>%
        dplyr::mutate(data_weight = prop_weight/dplyr::n())
    }
    
    # fit model 
    # - lmer v glmer - # - family
    if (is.null(fam)){
      if (incl_weights == F){
        tmp_m <- lme4::lmer(formula = orig_form, 
                            data = tr,
                            REML = REML)
      }
      if (incl_weights == T){
        tmp_m <- lme4::lmer(formula = orig_form, 
                            data = tr,
                            REML = REML,
                            weights = data_weight)
      }
    }
    if (!is.null(fam)){
      if (incl_weights == F){
        tmp_m <- lme4::glmer(formula = orig_form, 
                            data = tr,
                            family = fam)
      }
      if (incl_weights == T){
        tmp_m <- lme4::lmer(formula = orig_form, 
                            data = tr,
                            family = fam,
                            weights = data_weight)
      }
    }
    
    te$prediction <- predict(tmp_m, te, re.form = NA)
    
    cv_out_ls[[i]] <- te
    
    # get coefs, and cis...
    fixef_df <- as.data.frame(lme4::fixef(tmp_m))
    colnames(fixef_df) <- "Coefficient_value"
    fixef_df <- tibble::rownames_to_column(fixef_df, "Coefficient_name")
    # print(calc_CIs)
    if (!is.null(ci_type)){
      # message(crayon::green("Getting confidence intervals ..."))
      confint_df <- lme4::confint.merMod(tmp_m, parm = "beta_", method = ci_type)
      confint_df <- as.data.frame(confint_df)
      colnames(confint_df) <- c("loCI", "hiCI")
      confint_df$Coefficient_name <- rownames(confint_df)
      fixef_df <- dplyr::left_join(fixef_df,
                                   confint_df,
                                   by = "Coefficient_name")
    }
    
    fixef_df$fold_id <- fold_ids[i]
    coef_out_ls[[i]] <- fixef_df
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  cv_out_df <- dplyr::bind_rows(cv_out_ls)
  coef_out_df <- dplyr::bind_rows(coef_out_ls)
  
  # backtransform prediction, y as resp term, 
  if (!is.null(backtransform)){
    cv_out_df <- cv_out_df %>%
      dplyr::mutate(prediction = backtransform(prediction),
                    y = !!sym(all.vars(orig_form[[2]])))
  }
  # or convert resp value to y - based on eqn, and don't backtransform prediction
  if (is.null(backtransform)){
    cv_out_df <- cv_out_df %>%
      dplyr::mutate(y = !!orig_form[[2]] )
  }
  
  cv_out_df <- cv_out_df %>%
    dplyr::mutate(error = y - prediction)
  
  # TODO: 
  # cv metrics
  cv_metrics <- calc_cv_metrics(cv_out_df$prediction, 
                                cv_out_df$y, 
                                cv_out_df$error)
  
  return(list(cv_pred_df = cv_out_df, 
              cv_metrics = cv_metrics,
              cv_coef_df = coef_out_df,
              fold_type  = fold_type,
              fold_col   = fold_col,
              n_folds    = n_folds
              # ci_type    = ci_type
              ))
}

# hib_cv_inla

# hib_cv_brms


#' Calculate a range of metrics using cross-validation outputs
#' 
#' @description 
#' This function calculates a range of regression-focussed metrics to evaluate 
#' model performance on cross-validation output
#' @details 
#' Expects vectors of predictions, observations, and error
#' @param pred A [`numeric`] vector of predicted values obtained via cross-validation.
#' @param obs A [`numeric`] vector of observed values use in model-fitting
#' @param error A [`numeric`] vector of error (obs - pred)
#' @keywords model, evaluation
#' @examples 
#' \dontrun{
#'   cv_metrics <- calc_cv_metrics(predictions, observations, error)
#' } 
#' @noRd
calc_cv_metrics <- function(pred, obs, error){
  # Metrics - cor, cor^2, mse, rmse, bias (me), mae, predictive r2, rmsle (log pred and obs first), 
  # rse [relative sq error] (divide the MSE of model by the MSE of a model which uses 
  # the mean as the predicted value), rrmse (sqrt(mse/var))
  # nrmse (what base for normalisation??, sd)
  
  # checks - any nas...
  # align metrics with possible metrics..
  
  
  pred_mean <- mean(pred, na.rm = T)
  obs_mean <- mean(obs, na.rm = T)
  
  obs_minus_obs_mean <- obs - obs_mean
  
  # cor, corsq, 
  # Need to rm NAs...
  cor_out <- stats::cor(pred, obs)
  cor_sq <- cor_out^2
  
  # me, mae, rae, mape, mse, rmse, nrmse, rmsle, rse, rrmse
  me <- mean(error, na.rm = T)
  
  ae <- abs(error)
  
  mae <- mean(ae, na.rm = T)
  
  rae <- sum(ae, na.rm = T) / sum(abs(obs_minus_obs_mean), na.rm = T)
  
  # mape - issue of zeros in obs... - v small adjustment
  mape <- mean( (ae / obs), na.rm = T) * 100
  
  smape <- mean( (ae / ( (pred + obs) / 2) ), na.rm = T ) * 100
  
  mse <- mean(error^2, na.rm = T)
  rmse <- sqrt(mse)
  
  nrmse <- rmse / sd(pred, na.rm = T)
  
  # rmsle, rse, rrmse,
  rrmse <- sqrt( ( sum(error^2, na.rm = T) / sum(pred^2, na.rm = T) ) )
  
  rmsle <- sqrt(mean( ( log(obs + 1) - log(pred + 1) )^2 , na.rm = T))
  
  rse <- sum(error^2, na.rm = T) / sum( (obs - obs_mean)^2 , na.rm = T)
  
  # rsq - [1 - (PRESS / sums of squares total)] * 100
  rsq <- ( 1 - rse ) * 100
  
  # return list of vectors...
  metr <- c("Correlation" = cor_out,
            "Correlation_squared" = cor_sq,
            
            "Mean_error" = me, 
            "Mean_absolute_error" = mae, 
            "Relative_absolute_error" = rae, 
            "Mean_absolute_percentage_error" = mape, 
            "Symmetric_mean_absolute_percentage_error" = smape, 
            
            "Mean_squared_error" = mse, 
            "Root_mean_squared_error" = rmse, 
            "Normalised_root_mean_squared_error" = nrmse,
            "Relative_root_mean_squared_error" = rrmse,
            "Root_mean_squared_log_error" = rmsle,
            "Relative_squared_error" = rse, 
            
            "Predictive_R_squared" = rsq
            )
  
  return(metr)
}



