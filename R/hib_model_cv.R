# TODO: allow for custom N random - limit to 20? - have got 2,3,5,10 - most common


#' Cross-validate a fitted hib mixed-effects model
#' 
#' @description 
#' This function re-fits models on different subsets of the data, so that subsequent 
#' evaluation of model performance can take place.
#' @details 
#' Expects a hib model object...
#' @param mod_obj A [`hib_XX_model`] object.
#' @param db_obj A [`hib_database`] object, or, if modelling used multiple data-
#' sources, a list of database objects
#' @param fold_type A [`character`] indicating the type of cross-validation; one
#' of: "Random-2", "Random-3", "Random-5", "Random-10", "Source_ID", "SS", 
#'     "Country", "UN_subregion", "UN_region", "Biome", "Realm", "Custom"
#' @param custom_col A [`character`] indicating the column in the model data to 
#' use to define cross-validation folds, when fold_type == "Custom"
#' @param random_seed A [`numeric`] indicating the random seed used to define
#' random cross-validation folds - only used if 'Random-X' in fold_type
#' @param backtransform A function to convert the model prediction back to the original
#' scale of interest
#' @param use_original_mesh A logical value indicating whether cross-validation
#' of spatial INLA models should use the original mesh for all cv folds. Computationally
#' this is slightly more efficient, but may affect cv performance metrics - witheld
#' sites will still exist in the mesh. TODO: update - lme4_specs, brms_specs, INLA_specs
#' @keywords model, evaluation
#' @examples 
#' \dontrun{
#'   cv_out <- hib_mod_cv(mod_obj, db_obj, fold_type = "UN_region")
#' } 
#' @export
hib_mod_cv <- function(mod_obj, db_obj, fold_type, 
                       custom_col = NULL,
                       random_seed = NULL,
                       backtransform = NULL,
                       use_original_mesh = NULL){
  # TODO: option to just use original mesh for inla cv - override mesh_source
  
  # Assertion checks
  assertthat::assert_that(
    inherits(mod_obj, "hib_lme4_model") || inherits(mod_obj, "hib_inla_model") || inherits(mod_obj, "hib_brms_model"),     
    inherits(db_obj, "Database") || ( is.list(db_obj) & all(unlist(lapply(obj, function(x){inherits(x, "Database")}))) ),    
    is.character(fold_type),
    is.character(custom_col) || is.null(custom_col),
    is.numeric(random_seed) || is.null(random_seed),
    is.function(backtransform) || is.null(backtransform),
    is.logical(use_original_mesh) || is.null(use_original_mesh)
  )
  

  # Align input args with available options
  fold_type <- match.arg(fold_type,
                         c("Random-2", "Random-3", "Random-5", "Random-10", 
                           # "Random",
                           "Source_ID", "SS", #"SSBS",
                           "Country", "UN_subregion", "UN_region",
                           "Biome", "Realm",
                           "Custom"), # loo??
                         several.ok = T)
  fold_type <- unique(fold_type)
  
  # print(grepl("Random-", fold_type))
  # print(random_seed)
  if ((any(grepl("Random-", fold_type)) & is.null(random_seed)) ){
    random_seed <- 123
  }
  if ((!any(grepl("Random-", fold_type)) & !is.null(random_seed) )){
    random_seed <- NULL
    message(crayon::yellow(paste0("Warning: Random seed will not be used as ", 
                                  "'Random' cross-validation was not specified.")))
  }
  
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
  
  # make a spec df - each of fold types (unique), plus as many custom as unique
  # custom cols
  cv_spec <- data.frame(fold_type = fold_type[fold_type != "Custom"],
                        fold_col = fold_type[fold_type != "Custom"],
                        stringsAsFactors = F) |>
    dplyr::mutate(fold_col = dplyr::case_when(fold_type == "Random-2" ~ "fold_id",
                                              fold_type == "Random-3" ~ "fold_id",
                                              fold_type == "Random-5" ~ "fold_id",
                                              fold_type == "Random-10" ~ "fold_id",
                                              
                                              T ~ fold_col))
  
  
  if ("Custom" %in% fold_type){
    cv_spec <- bind_rows(cv_spec,
                         expand.grid(fold_type = c("Custom"),
                                     fold_col = custom_col,
                                     stringsAsFactors = F) |>
                           as.data.frame()
    )
  }
  

  if (inherits(mod_obj, "hib_inla_model")){
    # Package check
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop(
        "Package \"INLA\" must be installed to use this function.",
        call. = FALSE
      )
    }
    # OR
    # rlang::check_installed("INLA", reason = "to use INLA-based model fitting.")
    
    m_form <- as.character(mod_obj$user_formula$full)
    m_form <- paste(m_form[2], m_form[1], m_form[3])
    if (mod_obj$inla_incl_spatial_term){
      m_form <- paste(m_form,
                      "+ f(s, model = site.spde)")
    }
    # message(crayon::green("Cross-validating supplied model:\n\t" ))
    # message(crayon::green(format(formula(mod_obj$user_formula$full))))
    message(crayon::green(paste0("Cross-validating supplied model:\n\t", 
                                 m_form, "\n")))
  } else {
    message(crayon::green(paste0("Cross-validating supplied model:\n\t", 
                                 format(stats::formula(mod_obj$model)), "\n")))
  }
  
  
  if (inherits(mod_obj, "hib_brms_model")){
    message(crayon::green("Have some patience, this may take a while...\n"))
  }
  cv_out <- purrr::pmap(cv_spec,
                        function(fold_type, fold_col,  
                                 mod_ = mod_obj,
                                 db_ = db_obj,
                                 rand_seed = random_seed,
                                 backtr = backtransform,
                                 orig_mesh = use_original_mesh){
                          hib_cv_mod_inner(mod_, db_, fold_type, fold_col, 
                                           rand_seed, backtr, orig_mesh)
                        })
  
  # TODO: drop cv item if n_folds == 1
  # cv_out[[i]]n_folds
  cv_idx <- purrr::map(cv_out,
                       function(x){
                         x$n_folds > 1
                       }) |> 
    unlist()
  cv_out <- cv_out[cv_idx]
  
  message(crayon::green("Generating in-sample performance metrics for comparison.\n"))
  
  if (inherits(mod_obj, "hib_lme4_model")){
    # In-sample preds
    orig_form <- stats::formula(mod_obj$model)
    insample_df <- mod_obj$data
    insample_df$prediction <- predict(mod_obj$model,
                                      insample_df,
                                      re.form = NA)
  }
  
  # TODO:
  if (inherits(mod_obj, "hib_inla_model")){
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # message(crayon::yellow("Warning: INLA insample prediction in development, current output is preliminary!"))
    
    # In-sample preds
    insample_df <- mod_obj$data
    
    # Need orig_form (user supplied$response_fixed) to edit/backtransform response 
    orig_form <- mod_obj$user_formula$response_fixed
    fixed_covars <- all.vars(orig_form[3])
    
    fitted_form_ <- mod_obj$inla_fitted_formula
    orig_fam <- mod_obj$model$.args$family
    
    # Make a new inla stack for test predictions
    stk_te <- mod_obj$inla_stack
    
    # Drop y, rename index
    stk_te$data$data$y <- NA
    names(stk_te$data$index)[names(stk_te$data$index) == "fit"] <- "te"
    
    # Drop iid random effects
    rand_vars <- all.vars(mod_obj$user_formula$random)
    
    stk_te$effects$data <- stk_te$effects$data |> 
      dplyr::mutate_at(rand_vars, ~ NA)

    # print(head(stk_te$effects$data))
    # Combine stacks
    comb_stk <- INLA::inla.stack(mod_obj$inla_stack,
                                 stk_te)
    
    if (length(fixed_covars) > 0){
      check_flag <- 0
      while(check_flag == 0){
        tmp_m <- INLA::inla(formula = fitted_form_, # orig_form,
                            family = orig_fam,
                            data = INLA::inla.stack.data(comb_stk),
                            control.fixed = list(expand.factor.strategy='inla'),
                            # control.inla = list(strategy = "", int.strategy = ""), 
                            control.compute = list(dic = TRUE, waic = TRUE),
                            control.predictor = list(A = INLA::inla.stack.A(comb_stk)),
                            verbose = F)
        
        check_flag <- nrow(tmp_m$summary.fixed)
      }
    } else {
      tmp_m <- INLA::inla(formula = fitted_form_, # orig_form,
                          family = orig_fam,
                          data = INLA::inla.stack.data(comb_stk),
                          control.fixed = list(expand.factor.strategy='inla'),
                          # control.inla = list(strategy = "", int.strategy = ""), 
                          control.compute = list(dic = TRUE, waic = TRUE),
                          control.predictor = list(A = INLA::inla.stack.A(comb_stk)),
                          verbose = F)
    }
    
    
    # Extract prediction and cis to insample_df
    pred_idx <- INLA::inla.stack.index(comb_stk, tag = "te")$data
    tmp_pred <- tmp_m$summary.fitted.values[pred_idx,]
    
    insample_df$prediction <- tmp_pred$mean
    insample_df$prediction_q2.5 <- tmp_pred$`0.025quant`
    insample_df$prediction_q97.5 <- tmp_pred$`0.975quant`
    # Drop y - will be regenerated below for error estimation etc
    insample_df <- insample_df |> 
      dplyr::select(-y)
    # %%%%%
    # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  }
  
  if (inherits(mod_obj, "hib_brms_model")){
    # In-sample preds
    orig_form <- mod_obj$model$formula$formula
    insample_df <- mod_obj$data
    
    tmp_pred <- fitted(mod_obj$model,
                       insample_df,
                       re_form = NA) |>
      as.data.frame()
    insample_df$prediction <- tmp_pred$Estimate
    insample_df$prediction_q2.5 <- tmp_pred$Q2.5
    insample_df$prediction_q97.5 <- tmp_pred$Q97.5
    
  }
  
  # %%%%%
  # Insample processing - is the same for all model types..
  
  # backtransform prediction, y as resp term, 
  if (!is.null(backtransform)){
    insample_df <- insample_df |>
      dplyr::mutate(prediction = backtransform(prediction),
                    y = !!dplyr::sym(all.vars(orig_form[[2]])))
  }
  # or convert resp value to y - based on eqn, and don't backtransform prediction
  if (is.null(backtransform)){
    insample_df <- insample_df |>
      dplyr::mutate(y = !!orig_form[[2]] )
  }
  
  insample_df <- insample_df |>
    dplyr::mutate(error = y - prediction)
  
  insample_metrics <- calc_cv_metrics(insample_df$prediction, 
                                      insample_df$y, 
                                      insample_df$error)
  
  insample_res = list(insample_pred_df = insample_df, 
                      insample_metrics = insample_metrics)
  
  # %%%%%
  
  # %%%
  # Sort out use_original_mesh
  # if not an inla model, set to NA
  if (!inherits(mod_obj, "hib_inla_model")){
    use_original_mesh <- NA
  }

  # if an inla, and spatial model, check it's not null, if it is, NA, else logical T/F
  if (inherits(mod_obj, "hib_inla_model")){
    if (mod_obj$inla_incl_spatial_term == T){
      if (is.null(use_original_mesh)){
        use_original_mesh <- NA  
      }  
    }
    # if an inla, but non spatial model, set to na
    if (mod_obj$inla_incl_spatial_term == F){
      use_original_mesh <- NA
    }
  }
  # %%%
  
  # Return a cv object.. - incl original model? + df of cv predictions 
  # initialise cv_output...
  message(crayon::green("\nCollating cross-validation output ..."))
  
  cv_obj <- hib_cv_res$new(
    mod_obj = mod_obj,
    cv_res = cv_out,
    insamp_res = insample_res,
    use_orig_mesh = use_original_mesh)
  # Add use_original_mesh
  return(cv_obj)
}


#' Cross-validate a fitted hib_XX_model
#' 
#' @description 
#' This function re-fits the model on different subsets of the data, so that subsequent 
#' evaluation of model performance can take place.
#' @details 
#' Expects a hib_XX_model object...
#' @param mod_obj A [`hib_XX_model`] object.
#' @param db_obj A [`hib_database`] object, or, if modelling used multiple data-
#' sources, a list of database objects
#' @param fold_type A [`character`] indicating the type of cross-validation; one
#' of: ..., <spatial_options>, ..., "Custom" - ?"Random-5", "Random-10"?
#' @param fold_col A [`character`] indicating the column in the model data to 
#' use to define cross-validation folds
#' @param random_seed A [`numeric`] indicating the random seed used to define
#' random cross-validation folds - only used if 'Random-X' in fold_type
#' @param backtransform A function to convert the model prediction back to the original
#' scale of interest
#' @param use_original_mesh A logical value indicating whether cross-validation
#' of spatial INLA models should use the original mesh for all cv folds. Computationally
#' this is slightly more efficient, but may affect cv performance metrics - witheld
#' sites will still exist in the mesh.
#' @keywords model, evaluation
#' @examples 
#' \dontrun{
#'   cv_out <- hib_cv_mod_inner(mod_obj, db_obj, 
#'                          fold_type = "Source_ID",
#'                          fold_col  = "Source_ID")
#' } 
#' @noRd
hib_cv_mod_inner <- function(mod_obj, db_obj, 
                             fold_type, fold_col,
                             random_seed = NULL,
                             backtransform = NULL,
                             use_original_mesh = NULL){
  
  # Assertion checks 
  assertthat::assert_that(
    inherits(mod_obj, "hib_lme4_model") || inherits(mod_obj, "hib_inla_model") || inherits(mod_obj, "hib_brms_model"),
    inherits(db_obj, "Database") || is.list(db_obj), 
    is.character(fold_type),
    is.character(fold_col),
    is.numeric(random_seed) || is.null(random_seed),
    is.function(backtransform) || is.null(backtransform),
    is.logical(use_original_mesh) || is.null(use_original_mesh)
  )
  
  # Check fold_type matches one of the options
  fold_type <- match.arg(fold_type, 
                         c("Random-2", "Random-3", "Random-5", "Random-10", 
                           "Source_ID", "SS", #"SSBS",
                           "Country", "UN_subregion", "UN_region",
                           "Biome", "Realm",
                           "Custom"), # loo??
                         several.ok = F)
  
  # %%%%%
  fold_type_msg <- fold_type
  if (fold_type == "Custom"){
    fold_type_msg <- paste0(fold_type, " (", fold_col, ")")
  }
  message(crayon::green(paste0("\nRunning cross-validation using: ", fold_type_msg, "\n")))
  
  
  # Ensure fold_col is in database, and extract ...
  if (fold_type %in% c("Source_ID", "SS", "Country", "UN_subregion", "UN_region",
                       "Biome", "Realm", "Custom")){
    # Check fold_col is in db_obj$site
    if(inherits(db_obj, "Database")){
      assertthat::assert_that(
        assertthat::has_name(db_obj$sites, fold_col)
      )
      fold_df <- db_obj$sites[,c("SSBS", fold_col)]
      
      # TODO: what if model is of beta diversity - need to use SSBS1? 
      # - here, assuming that SSBS1/2 are a pair, with shuffling based on SSBS1
      # TODO: get lat/long of each ssbs
    }
    # or each db_obj$site if db_obj is a list
    if(is.list(db_obj)){
      assertthat::assert_that(
        all(unlist(lapply(db_obj, function(x){fold_col %in% colnames(x$sites)})))
      )
      
      fold_df <- lapply(db_obj, function(x){x$sites[,c("SSBS", fold_col)]}) |>
        dplyr::bind_rows()
      
      # TODO: what if model is of beta diversity - need to use SSBS1?
      # TODO: get lat/long of each ssbs (midpoint of SS?)
    }
  }
  
  
  # %%%%%%%%%%
  # Random fold generation
  if (fold_type %in% c("Random-2", "Random-3", "Random-5", "Random-10" )){
    if ("SSBS" %in% colnames(mod_obj$data)){
      set.seed(random_seed)
      fold_df <- data.frame(SSBS = sample(mod_obj$data$SSBS))  
    }
    if ("SSBS1" %in% colnames(mod_obj$data)){
      set.seed(random_seed)
      fold_df <- data.frame(SSBS = sample(mod_obj$data$SSBS1))  
      # TODO: is this a good idea? - maybe should sample from SSBS1_SSBS2 - to 
      # properly "split" data
      # Especially as SSBS1 will be duplicated many times in data - mostly going
      # to be primary sites - this random sampling could be messy
    }
    
    n_rows <- nrow(mod_obj$data)
    if (fold_type == "Random-2"){
      fold_id_vec <- get_rand_fold_id(n_rows, 2)
    }
    if (fold_type == "Random-3"){
      fold_id_vec <- get_rand_fold_id(n_rows, 3)
    }
    if (fold_type == "Random-5"){
      fold_id_vec <- get_rand_fold_id(n_rows, 5)
    }
    if (fold_type == "Random-10"){
      fold_id_vec <- get_rand_fold_id(n_rows, 10)
    }
    fold_df$fold_id <- fold_id_vec
  }
  # %%%%%%%%%%
  
  
  
  # NOTE: this data prep stuff should only happen once, before any cv_inner,
  # info is relevant to all cv types!!
  
  # Should we also try to "keep" SS together in other blockings...
  
  mod <- mod_obj$model
  orig_data <- mod_obj$data
  
  incl_weights <- F
  if ("data_weight" %in% colnames(orig_data)){
    incl_weights <- T
  }
  
  if (inherits(mod_obj, "hib_lme4_model")){
    # %%%%%%%%%%
    # LME4-specific data prep 
    # %%%%%%%%%%
    orig_form <- stats::formula(mod)
    
    fam <- NULL
    if (is(mod@resp, "glmResp")){
      fam <- mod@resp$family
    }
    
    if (mod@devcomp$dims[["REML"]]){
      REML <- T
    } else {REML <- F}
    
    ci_type <- mod_obj$ci_type
    # %%%%%%%%%%
    # LME4-specific dataprep - END
    # %%%%%%%%%%
  }
  
  if (inherits(mod_obj, "hib_inla_model")){
    # %%%%%%%%%%
    # INLA-specific dataprep 
    # %%%%%%%%%%
    
    # Formulas 
    
    orig_form <- mod_obj$user_formula$full
    fitted_form <- as.formula(mod$.args$formula)
    fitted_form_ <- mod_obj$inla_fitted_formula
    
    orig_form_user_ls <- mod_obj$user_formula
    # $full
    # $response_fixed
    # $random
    
    # family
    orig_fam <- mod$.args$family
    
    # Spatial settings
    orig_incl_spatial_term <- NA
    mesh_source <- NA
    orig_mesh <- NA
    orig_mesh_settings <- NA
    
    orig_incl_spatial_term <- mod_obj$inla_incl_spatial_term
    if (orig_incl_spatial_term){
      if (is.null(use_original_mesh)){
        use_original_mesh <- NA
      }
      # is mesh user provided or based on mesh_settings?
      mesh_source <- mod_obj$inla_mesh_source
      if (mesh_source == "user_provided"){
        orig_mesh <-  mod_obj$inla_mesh
      }
      if (mesh_source == "mesh_settings"){
        orig_mesh_settings <-  mod_obj$inla_mesh_settings
      }
      if (use_original_mesh == T){
        orig_mesh <-  mod_obj$inla_mesh
        mesh_source <- "user_provided"
      }
    }
    
    # prior_list ...
    # fit_settings ...
    # %%%%%%%%%%
    # INLA-specific dataprep -  END
    # %%%%%%%%%%
  }
  
  if (inherits(mod_obj, "hib_brms_model")){
    # %%%%%%%%%%
    # BRMS-specific dataprep -  
    # %%%%%%%%%%
    orig_form <- mod$formula$formula
    orig_fam <- mod$family
    orig_extra_args <- mod_obj$extra_args
    
    m_args <- orig_extra_args
    m_args[["formula"]] <- orig_form
    m_args[["family"]] <- orig_fam
    # %%%%%%%%%%
    # BRMS-specific dataprep -  END
    # %%%%%%%%%%
  }
  
  
  # %%%%%%%%%%%%%%%%%%
  # Below is specific to each cv type... - i.e. needed per cv_inner
  
  
  # If model is based on contrasts (i.e., comp-sim style), need to get, info on both
  # SSBS1 and SSBS2... - below, am just using SSBS1...
  
  # Add important fold variable - from db_obj - based on SSBS, or, 
  # if SSBS1/SSBS2 in orig_data, use these...
  if ("SSBS" %in% colnames(orig_data)){
    # TODO: BUG FIX NEEDED: 20231005 - if fold column is already in orig_data, 
    # Error occurs as you then get fold_col.x and fold_col.y
    # If only col in fold_df (other than SSBS/SSBS1) is already present, 
    # drop and re add, or just don't add? - question re SSBS1... - ether seems okay
    # If already present, just don't add, assume model df has correct var...
    if (fold_col %in% colnames(orig_data)){
      cv_data <- orig_data
    }
    if (!fold_col %in% colnames(orig_data)){
      cv_data <- dplyr::left_join(orig_data,
                                  fold_df,
                                  by = "SSBS")
    }
    # cv_data <- dplyr::left_join(orig_data,
    #                             fold_df,
    #                             by = "SSBS")
  }
  if ("SSBS1" %in% colnames(orig_data)){
    if (fold_col %in% colnames(orig_data)){
      cv_data <- orig_data
    }
    if (!fold_col %in% colnames(orig_data)){
      cv_data <- dplyr::left_join(orig_data,
                                  fold_df,
                                  by = c("SSBS1" = "SSBS"))
    }
  }
  
  fold_ids <- unique(cv_data |> dplyr::pull(!!fold_col))
  n_folds <- length(fold_ids)
  # print(n_folds)
  # Make a list to store output - length = length(unique(new_df[,fold_var]))
  cv_out_ls <- vector("list", n_folds)
  coef_out_ls <- vector("list", n_folds)
  
  # TODO: if only one fold !! - need to note this...
  if (n_folds == 1){
    message(crayon::yellow(paste0("Warning: ", fold_type_msg, " has only 1 fold... Cross-validation for this specification will not be run.\n")))
    
    return(list(cv_pred_df  = NA, 
                cv_metrics  = NA,
                cv_fold_metrics  = NA,
                cv_coef_df  = NA,
                fold_type   = fold_type,
                fold_col    = fold_col,
                n_folds     = n_folds,
                random_seed = NA
    ))
  }
  
  message(crayon::green(paste0("\tTotal number of folds: ", n_folds, "\n\n\tProgress")))
  
  # TODO: this is only really good for LME4 models, unless we make brms/inla
  # silent
  if ( (inherits(mod_obj, "hib_lme4_model")) | 
        (inherits(mod_obj, "hib_inla_model")) ){
    # print("WOOP - lme4/inla")
    pb <- txtProgressBar(min = 1, max = max(2, n_folds), style = 3)
  }
  
  for (i in 1:n_folds){
    
    tr <- cv_data |>
      dplyr::filter(!!dplyr::sym(fold_col) != fold_ids[i])
    te <- cv_data |>
      dplyr::filter(!!dplyr::sym(fold_col) == fold_ids[i])
    
    
    # Make a fold_id/fold_idx column for te...
    te <- te |>
      dplyr::mutate(fold_id = !!dplyr::sym(fold_col))
    
    # Adjust weights in tr so sum per data type - dataset_modelling_id
    # sum to same value as in orig data...
    if (incl_weights == T){
      weight_df <- orig_data |>
        dplyr::group_by(dataset_modelling_id) |>
        dplyr::summarise(prop_weight = sum(data_weight))
      
      tr <- tr |>
        dplyr::select(-data_weight)
      te <- te |>
        dplyr::select(-data_weight)
      
      tr <- dplyr::left_join(tr, weight_df) |>
        dplyr::group_by(dataset_modelling_id) |>
        dplyr::mutate(data_weight = prop_weight/dplyr::n())
    }
    
    if (inherits(mod_obj, "hib_lme4_model")){
      # %%%%%%%%%%
      # LME4-specific
      # %%%%%%%%%%
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
          tmp_m <- lme4::glmer(formula = orig_form, 
                               data = tr,
                               family = fam,
                               weights = data_weight)
        }
      }
      
      te$prediction <- predict(tmp_m, te, re.form = NA)
      
      # get coefs, and cis...
      fixef_df <- as.data.frame(lme4::fixef(tmp_m))
      colnames(fixef_df) <- "Coefficient_value"
      fixef_df <- tibble::rownames_to_column(fixef_df, "Coefficient_name")
      
      if (!is.null(ci_type)){
        confint_df <- lme4::confint.merMod(tmp_m, parm = "beta_", method = ci_type)
        confint_df <- as.data.frame(confint_df)
        colnames(confint_df) <- c("loCI", "hiCI")
        confint_df$Coefficient_name <- rownames(confint_df)
        fixef_df <- dplyr::left_join(fixef_df,
                                     confint_df,
                                     by = "Coefficient_name")
      }
      # %%%%%%%%%%
      # LME4-specific - END
      # %%%%%%%%%%
    }
    
    if (inherits(mod_obj, "hib_inla_model")){
      # %%%%%%%%%%
      # INLA-specific
      # %%%%%%%%%%
      
      # Make X_tr and X_te
      X_tr <- tr |> 
        dplyr::select(-y)
      
      # Drop iid random effects from te
      rand_vars <- all.vars(orig_form_user_ls$random)
    
      X_te <- te |> 
        dplyr::select(-y) |>
        dplyr::mutate_at(rand_vars, ~ NA)
      
      # print(dim(X_tr))
      # print(dim(X_te))
      
      # If not a spatial model
      if (orig_incl_spatial_term == F){
        # Just use model dfs
        stk_tr <- INLA::inla.stack(tag = "tr",
                                   data = list(y = tr$y),
                                   A = list(1),
                                   effects = list(X = X_tr))
        
        stk_te <- INLA::inla.stack(tag = "te",
                                   data = list(y = NA),
                                   A = list(1),
                                   effects = list(X = X_te))
      }
      
      # TODO:
      # If spatial, make and generate mesh...
      if (orig_incl_spatial_term == T){
        
        # print(orig_mesh_settings)
        
        # Get lat/long
        if (("Longitude" %in% colnames(X_tr)) &
            ("Latitude" %in% colnames(X_tr)) ){
          coords_tr <- as.matrix(X_tr[,c("Longitude", "Latitude")])
          coords_te <- as.matrix(X_te[,c("Longitude", "Latitude")])
        }
        if (("Longitude_midpoint" %in% colnames(X_tr)) &
            ("Latitude_midpoint" %in% colnames(X_tr)) ){
          coords_tr <- as.matrix(X_tr[,c("Longitude_midpoint", "Latitude_midpoint")])
          coords_te <- as.matrix(X_te[,c("Longitude_midpoint", "Latitude_midpoint")])
        }
      
        
        if (mesh_source == "mesh_settings"){
          # Need to re-build a mesh - based only on tr points
          # Spatial settings
          # orig_mesh
          # orig_mesh_settings
          
          if (orig_mesh_settings$mesh_dim == "3d"){
            # print("Building 3d mesh...")
            coords_tr <- INLA::inla.mesh.map(coords_tr, projection = "longlat")
            cv_mesh <- INLA::inla.mesh.create(loc = coords_tr,
                                                 globe = 8)
          }
          
          if (orig_mesh_settings$mesh_dim == "2d"){
            # print("Building 2d mesh...")
            # TODO: Why is building the cv mesh taking soooo long...
            cv_mesh <- INLA::inla.mesh.2d(loc = coords_tr,
                                             max.edge = orig_mesh_settings$max.edge,
                                             cutoff = orig_mesh_settings$cutoff,
                                             offset = orig_mesh_settings$offset)
            # print("Mesh built!")
          }
        }
        
        if (mesh_source == "user_provided"){
          cv_mesh <- orig_mesh
          # TODO: CHECK THIS IS WORKING
        }
        
        
        # spde
        # Set up spde model
        site.spde <- INLA::inla.spde2.matern(cv_mesh, alpha = 2)
        s.site <- INLA::inla.spde.make.index('s', n.spde = site.spde$n.spde)
        # TODO: add prior specs...
        # TODO: pcmatern ??
        
        # Create weights (A) matrix
        # -for tr
        A_tr <- INLA::inla.spde.make.A(cv_mesh, 
                                       loc = coords_tr)
        # -for te
        A_te <- INLA::inla.spde.make.A(cv_mesh, 
                                       loc = coords_te)
        
        # Combine into stacks - tr and te
        stk_tr <- INLA::inla.stack(tag = "tr",
                                   data = list(y = tr$y),
                                   A = list(1, A_tr),
                                   effects = list(X = X_tr,
                                                  s = s.site))
        
        stk_te <- INLA::inla.stack(tag = "te",
                                   data = list(y = NA),
                                   A = list(1, A_te),
                                   effects = list(X = X_te,
                                                  s = s.site)) 
      }
      
      # Combine stacks
      comb_stk <- INLA::inla.stack(stk_tr, 
                                   stk_te)
      # Fit model
      fixed_covars <- all.vars(orig_form_user_ls$response_fixed[[3]])
      
      if (length(fixed_covars) > 0){
        check_flag <- 0
        while(check_flag == 0){
          tmp_m <- INLA::inla(formula = fitted_form_, # orig_form,
                           family = orig_fam,
                           data = INLA::inla.stack.data(comb_stk),
                           control.fixed = list(expand.factor.strategy='inla'),
                           # control.inla = list(strategy = "", int.strategy = ""), 
                           control.compute = list(dic = TRUE, waic = TRUE),
                           control.predictor = list(A = INLA::inla.stack.A(comb_stk)),
                           verbose = F)
          
          check_flag <- nrow(tmp_m$summary.fixed)
        }
      } else {
        tmp_m <- INLA::inla(formula = fitted_form_, # orig_form,
                         family = orig_fam,
                         data = INLA::inla.stack.data(comb_stk),
                         control.fixed = list(expand.factor.strategy='inla'),
                         # control.inla = list(strategy = "", int.strategy = ""), 
                         control.compute = list(dic = TRUE, waic = TRUE),
                         control.predictor = list(A = INLA::inla.stack.A(comb_stk)),
                         verbose = F)
      }
      
      
      # Extract prediction and cis to te (a df of data + pred)
      pred_idx <- INLA::inla.stack.index(comb_stk, tag = "te")$data
      tmp_pred <- tmp_m$summary.fitted.values[pred_idx,]
      
      te$prediction <- tmp_pred$mean
      te$prediction_q2.5 <- tmp_pred$`0.025quant`
      te$prediction_q97.5 <- tmp_pred$`0.975quant`
      # Drop y from te - will be regenerated below for error estimation etc
      te <- te |> 
        dplyr::select(-y)
      
      # Get coefficients from model
      fixef_df <- tmp_m$summary.fixed |>
        tibble::rownames_to_column("Coefficient_name") |>
        dplyr::rename(Coefficient_value = mean,
                      loCI = `0.025quant`,
                      median = `0.5quant`,
                      hiCI = `0.975quant`)
      
      # %%%%%%%%%%
      # INLA-specific - END
      # %%%%%%%%%%
    }
    
    if (inherits(mod_obj, "hib_brms_model")){
      # %%%%%%%%%%
      # BRMS-specific
      # %%%%%%%%%%
      
      # fit model 
      tmp_args <- NULL
      tmp_args <- m_args
      tmp_args[["data"]] <- tr
      
      tmp_m <- do.call(brms::brm,
                       tmp_args)
      
      # predict or fitted??
      # tmp_pred <- predict(tmp_m, te, re_form = NA)|>
      #   as.data.frame()
      tmp_pred <- fitted(tmp_m, te, re_form = NA) |>
        as.data.frame()
      te$prediction <- tmp_pred$Estimate
      te$prediction_q2.5 <- tmp_pred$Q2.5
      te$prediction_q97.5 <- tmp_pred$Q97.5
      
      # get coefs, and cis...
      fixef_df <- as.data.frame(brms::fixef(tmp_m))
      colnames(fixef_df) <- c("Coefficient_value", "Coefficient_error", "loCI", "hiCI")
      fixef_df <- tibble::rownames_to_column(fixef_df, "Coefficient_name")
      
      # %%%%%%%%%%
      # BRMS-specific - END
      # %%%%%%%%%%
    }
    
    cv_out_ls[[i]] <- te
    
    fixef_df$fold_id <- fold_ids[i]
    coef_out_ls[[i]] <- fixef_df
    
    if ( (inherits(mod_obj, "hib_lme4_model")) |
         (inherits(mod_obj, "hib_inla_model")) ){
      setTxtProgressBar(pb, i)
    }
  }
  
  if ( (inherits(mod_obj, "hib_lme4_model")) |
        (inherits(mod_obj, "hib_inla_model")) ){
    close(pb)
  }
  
  cv_out_df <- dplyr::bind_rows(cv_out_ls)
  coef_out_df <- dplyr::bind_rows(coef_out_ls)
  
  # backtransform prediction, y as resp term, 
  if (!is.null(backtransform)){
    cv_out_df <- cv_out_df |>
      dplyr::mutate(prediction = backtransform(prediction),
                    y = !!dplyr::sym(all.vars(orig_form[[2]])))
  }
  # or convert resp value to y - based on eqn, and don't backtransform prediction
  if (is.null(backtransform)){
    cv_out_df <- cv_out_df |>
      dplyr::mutate(y = !!orig_form[[2]] )
  }
  
  cv_out_df <- cv_out_df |>
    dplyr::mutate(error = y - prediction)
  
  # cv metrics
  cv_metrics <- calc_cv_metrics(cv_out_df$prediction, 
                                cv_out_df$y, 
                                cv_out_df$error)
  
  cv_fold_metrics <- calc_per_fold_cv_metrics(cv_out_df$prediction, 
                                              cv_out_df$y, 
                                              cv_out_df$error,
                                              cv_out_df$fold_id)
  
  # Only return random seed if random-
  if (fold_type %in% c(c("Random-2", "Random-3", "Random-5", "Random-10" ))){
    rand_num <- random_seed
  } else {
    rand_num <- NA
  }
  
  return(list(cv_pred_df  = cv_out_df, 
              cv_metrics  = cv_metrics,
              cv_fold_metrics  = cv_fold_metrics,
              cv_coef_df  = coef_out_df,
              fold_type   = fold_type,
              fold_col    = fold_col,
              n_folds     = n_folds,
              random_seed = rand_num
  ))
}


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
  cor_out <- stats::cor(pred, obs, use = "complete.obs")
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
  
  rrmse <- sqrt( ( sum(error^2, na.rm = T) / sum(pred^2, na.rm = T) ) )
  
  rmsle <- sqrt(mean( ( log(obs + 1) - log(pred + 1) )^2 , na.rm = T))
  
  rse <- sum(error^2, na.rm = T) / sum( obs_minus_obs_mean^2 , na.rm = T)
  
  # rsq - [1 - (PRESS / sums of squares total)] * 100
  rsq <- ( 1 - rse )
  
  # TODO:
  # Any metric that is relative to "obs_mean" or similar could also be calculated per fold 
  # - global mean is not per-fold training sample mean...
  # rsq, rse, rae - all relative to using mean as a predictor...
  # need fold_id...
  
  # CIs for mse, rmse, me...- sheiner and beale 1981
  n_ <- sum(!is.na(error))
  se_me <- sqrt( (1 / (n_ * (n_ - 1) ) ) * sum( (error - me)^2, na.rm = T) )
  se_mse <- sqrt( (1 / (n_ * (n_ - 1) ) ) * sum( (error^2 - mse)^2, na.rm = T) )
  
  me_lo <- me - (1.96 * se_me)
  me_hi <- me + (1.96 * se_me)
  
  mse_lo <- mse - (1.96 * se_mse)
  mse_hi <- mse + (1.96 * se_mse)
  
  if (mse_lo < 0){
    mse_lo <- 0
  }

  rmse_lo <- sqrt(mse_lo)
  rmse_hi <- sqrt(mse_hi)
  
  # Make metr into a df
  metr_df <- data.frame(Metric = c("Correlation",
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
                                   "Predictive_R_squared"),
                        Metric_val = c(cor_out, cor_sq,
                                       me, mae, rae, mape, smape, 
                                       mse, rmse, nrmse, 
                                       rrmse, rmsle, rse, 
                                       rsq),
                        Metric_lo = c(NA, NA, 
                                      me_lo, NA, NA, NA, NA,
                                      mse_lo, rmse_lo, NA,
                                      NA, NA, NA, 
                                      NA),
                        Metric_hi = c(NA, NA, 
                                      me_hi, NA, NA, NA, NA,
                                      mse_hi, rmse_hi, NA,
                                      NA, NA, NA, 
                                      NA),
                        Metric_shrt = c("cor", "corsq",
                                        "me", "mae", "rae", "mape", 
                                        "smape", 
                                        "mse", "rmse", "nrmse", 
                                        "rrmse", "rmsle", "rse", 
                                        "rsq")) |>
    dplyr::mutate(Metric_fct = factor(Metric,
                                      levels = c("Correlation",
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
                                                 "Predictive_R_squared"),
                                      labels = c("Correlation",
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
                                                 "Predictive~R^2")))
  
  
  return(metr_df)
}


# %%%%
#' Calculate a range of per-fold metrics using cross-validation outputs
#' 
#' @description 
#' This function calculates a range of  metrics to evaluate 
#' model performance on cross-validation output
#' @details 
#' Expects vectors of predictions, observations, error, and fold_id
#' @param pred A [`numeric`] vector of predicted values obtained via cross-validation.
#' @param obs A [`numeric`] vector of observed values use in model-fitting
#' @param error A [`numeric`] vector of error (obs - pred)
#' @param fold_id A [`numeric`] vector of error (obs - pred)
#' @keywords model, evaluation
#' @examples 
#' \dontrun{
#'   cv_metrics <- calc_cv_metrics(predictions, observations, error)
#' } 
#' @noRd
calc_per_fold_cv_metrics <- function(pred, obs, error, fold_id){
  # Metrics - cor, cor^2, mse, rmse, bias (me), mae, predictive r2, rmsle (log pred and obs first), 
  # rse [relative sq error] (divide the MSE of model by the MSE of a model which uses 
  # the mean as the predicted value), rrmse (sqrt(mse/var))
  # nrmse (what base for normalisation??, sd)
  
  # checks - any nas...
  # align metrics with possible metrics..
  
  # rae, rse, rsq
  tmp <- data.frame(pred    = pred, 
                    obs     = obs, 
                    error   = error, 
                    fold_id = fold_id)
  
  metr_df <- tmp |>
    dplyr::group_by(fold_id) |>
    dplyr::summarise(pred_mean = mean(pred, na.rm = T),
                     obs_mean  = mean(obs, na.rm = T),
                     
                     # %%%
                     # cor, corsq, 
                     # Need to rm NAs...
                     cor_out = stats::cor(pred, obs, use = "complete.obs"),
                     cor_sq  = cor_out^2,
                     # me, mae, rae, mape, mse, rmse, nrmse, rmsle, rse, rrmse
                     me = mean(error, na.rm = T),
                     mae = mean(abs(error), na.rm = T),
                     # mape - issue of zeros in obs... - v small adjustment
                     mape = mean( (abs(error) / obs), na.rm = T) * 100,
                     smape = mean( (abs(error) / ( (pred + obs) / 2) ), na.rm = T ) * 100,
                     mse = mean(error^2, na.rm = T),
                     rmse = sqrt(mse),
                     nrmse = rmse / sd(pred, na.rm = T),
                     rrmse = sqrt( ( sum(error^2, na.rm = T) / sum(pred^2, na.rm = T) ) ),
                     rmsle = sqrt(mean( ( log(obs + 1) - log(pred + 1) )^2 , na.rm = T)),
                     # %%%
                     
                     rae = sum(abs(error), na.rm = T) / sum((obs - obs_mean), na.rm = T),
                     rse = sum(error^2, na.rm = T) / sum( (obs - obs_mean)^2 , na.rm = T),
                     rsq = 1- rse 
                     ) |>
    dplyr::select(fold_id, cor_out, cor_sq, me, mae, mape, smape, mse, rmse, nrmse,
                  rrmse, rmsle,
                  rae, rse, rsq) |>
    tidyr::pivot_longer(!fold_id, 
                        names_to = "Metric_shrt", 
                        values_to = "Metric_val") |>
    dplyr::mutate(Metric_shrt = dplyr::case_when(Metric_shrt == "cor_out" ~ "cor",
                                                  Metric_shrt == "cor_sq" ~ "corsq",
                                                  T ~ Metric_shrt)) |>
    dplyr::mutate(Metric = dplyr::case_when(Metric_shrt == "cor"   ~ "Correlation",
                                            Metric_shrt == "corsq" ~ "Correlation_squared",
                                            Metric_shrt == "me"    ~ "Mean_error", 
                                            Metric_shrt == "mae"   ~ "Mean_absolute_error", 
                                            Metric_shrt == "rae"   ~ "Relative_absolute_error", 
                                            Metric_shrt == "mape"  ~ "Mean_absolute_percentage_error", 
                                            Metric_shrt == "smape" ~ "Symmetric_mean_absolute_percentage_error", 
                                            Metric_shrt == "mse"   ~ "Mean_squared_error", 
                                            Metric_shrt == "rmse"  ~ "Root_mean_squared_error", 
                                            Metric_shrt == "nrmse" ~ "Normalised_root_mean_squared_error",
                                            Metric_shrt == "rrmse" ~ "Relative_root_mean_squared_error",
                                            Metric_shrt == "rmsle" ~ "Root_mean_squared_log_error",
                                            Metric_shrt == "rse"   ~ "Relative_squared_error", 
                                            Metric_shrt == "rsq"   ~ "Predictive_R_squared",
                                            T ~ NA),
                  Metric_fct = factor(Metric,
                                      levels = c("Correlation",
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
                                                 "Predictive_R_squared"),
                                      labels = c("Correlation",
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
                                                 "Predictive~R^2")))
  
  return(metr_df)
  
}
# %%%%


#' Helper function to generate roughly evenly sized random folds
#' 
#' @description 
#' This function returns a vector of length n_obs with the values 1:n_folds
#' @details 
#' Expects numeric inputs
#' @param n_obs A [`numeric`] value indicating the number of data points
#' @param n_folds A [`numeric`] value indicating the number of folds required
#' @keywords model, evaluation
#' @examples 
#' \dontrun{
#'   fold_ids <- get_rand_fold_id(nrow(df), 5)
#' } 
#' @noRd
get_rand_fold_id <- function(n_obs, n_folds){
  assertthat::assert_that(
    is.numeric(n_obs),
    length(n_obs) == 1,
    is.numeric(n_folds),
    length(n_folds) == 1
  )
  if (n_obs %% n_folds == 0){
    fold_id_vec <- rep(1:n_folds, n_obs / n_folds )
  } else {
    fold_id_vec <- c( rep(1:n_folds, floor(n_obs / n_folds) ), 1:(n_obs %% n_folds) )
  }
  
  assertthat::assert_that(
    length(fold_id_vec) == n_obs
  )
  
  return(fold_id_vec)
}