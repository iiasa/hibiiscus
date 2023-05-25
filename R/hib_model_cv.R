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
#' of: ..., <spatial_options>, ..., "custom" - ?"random-5", "random-10"?
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
    inherits(db_obj, "Database") || is.list(db_obj),
    is.character(fold_type),
    is.character(custom_col) || is.null(custom_col),
    is.function(backtransform) || is.null(backtransform)
  )
  if (inherits(mod_obj, "hib_lme4_model")){
    cv_out <- hib_cv_lme4(mod_obj, db_obj, fold_type, custom_col,
                          backtransform)
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
  # key output cols - prediction, error, fold_id, fold_idx
  return(cv_out)
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
#' of: ..., <spatial_options>, ..., "custom" - ?"random-5", "random-10"?
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
hib_cv_lme4 <- function(mod_obj, db_obj, fold_type, custom_col = NULL,
                        backtransform = NULL){
  
  # Assertion checks 
  assertthat::assert_that(
    inherits(mod_obj, "hib_lme4_model"),
    inherits(db_obj, "Database") || is.list(db_obj), # TODO list of db...
    is.character(fold_type),
    is.character(custom_col) || is.null(custom_col),
    is.function(backtransform) || is.null(backtransform)
  )
  
  
  # Check fold_type matches one of the options
  fold_type <- match.arg(fold_type, 
                         c("random-5", "random-10", 
                           "Source_ID", "SS",
                           "Country", "UN_subregion", "UN_region",
                           "Biome", "Realm",
                           "custom"), # loo??
                         several.ok = F)
  fold_col <- NA
  if (fold_type %in% c("Source_ID", "SS", "Country", "UN_subregion", "UN_region",
                       "Biome", "Realm")){
    fold_col <- fold_type
  }
  if (fold_type == "custom"){
    fold_col <- custom_col
  }
  
  
  # Ensure fold_col is in database, and extract ...
  if (fold_type %in% c("Source_ID", "SS", "Country", "UN_subregion", "UN_region",
                       "Biome", "Realm", "custom")){
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
  if (fold_type == "random-5"){
    # Randomly sample/shuffle SSBS, assign to 5 indices
    fold_df <- data.frame(SSBS = sample(mod_obj$data$SSBS),
                          fold_id = seq(1,5))
    fold_col <- "fold_id"
  }

  if (fold_type == "random-10"){
    # Randomly sample/shuffle SSBS, assign to 10 indices
    fold_df <- data.frame(SSBS = sample(mod_obj$data$SSBS),
                          fold_id = seq(1,10))
    fold_col <- "fold_id"
  }
  
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
                                fold_df)
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
    
  # For every fold - subset data
  for (i in 1:n_folds){
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
      message(crayon::green("Getting confidence intervals ..."))
      confint_df <- lme4::confint.merMod(tmp_m, parm = "beta_", method = ci_type)
      confint_df <- as.data.frame(confint_df)
      colnames(confint_df) <- c("loCI", "hiCI")
      confint_df$Coefficient_name <- rownames(confint_df)
      fixef_df <- dplyr::left_join(fixef_df,
                                   confint_df)
    }
    
    fixef_df$fold_id <- fold_ids[i]
    coef_out_ls[[i]] <- fixef_df
  }
  
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
  
  return(list(cv_pred_df = cv_out_df, 
              cv_coef_df = coef_out_df))
}

# hib_cv_inla

# hib_cv_brms




