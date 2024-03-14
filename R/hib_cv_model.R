
#' Cross-validate a fitted hib mixed-effects model
#' 
#' @description 
#' This function re-fits models on different subsets of the data, so that subsequent 
#' evaluation of model performance can take place.
#' @details 
#' Expects a HibModel.
#' @param mod_obj One of; [`HibModelLME4`], [`HibModelINLA`] or [`HibModelBRMS`].
#' @param db_obj A [`HibDatabase`] object, or, if modelling used multiple data-
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
#' @param re_form A [`list`] indicating the random effects variables to be
#' included in the cross-validation/insample predicitons - if multiple items supplied,
#' each option will be used in a separate cv prediction.
#' Values can be either NA (no random effects), NULL (all random effects), or a
#' rhs formula specifying particular random effects to include.
#' Note: only test set data points with random-effect values present in the 
#' training data will have predictions if re_form is not NA.
#' @param use_original_mesh A logical value indicating whether cross-validation
#' of spatial INLA models should use the original mesh for all cv folds. Computationally
#' this is slightly more efficient, but may affect cv performance metrics - witheld
#' sites will still exist in the mesh. 
#' @keywords model, evaluation
#' @examples 
#' \dontrun{
#'   cv_out <- hib_mod_cv(mod_obj, db_obj, fold_type = "UN_region")
#' } 
#' @export
hib_cv_model <- function(mod_obj, db_obj, fold_type, 
                       custom_col = NULL,
                       random_seed = NULL,
                       re_form = list(NA),
                       backtransform = NULL,
                       use_original_mesh = NULL){
  # hib_mod_cv
  # TODO: option to just use original mesh for inla cv - override mesh_source
  # TODO: custom_error_func = list(NA)
  # should be a named list list("Metric_name" = function)

  # Assertion checks
  assertthat::assert_that(
    inherits(mod_obj, "HibModelLME4") || inherits(mod_obj, "HibModelINLA") || inherits(mod_obj, "HibModelBRMS"),     
    inherits(db_obj, "HibDatabase") || ( is.list(db_obj) & all(unlist(lapply(db_obj, function(x){inherits(x, "HibDatabase")}))) ),    
    is.character(fold_type),
    is.character(custom_col) || is.null(custom_col),
    is.numeric(random_seed) || is.null(random_seed),
    is.list(re_form), 
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
    cv_spec <- dplyr::bind_rows(cv_spec,
                         expand.grid(fold_type = c("Custom"),
                                     fold_col = custom_col,
                                     stringsAsFactors = F) |>
                           as.data.frame()
    )
  }
  
  if (inherits(mod_obj, "HibModelLME4")){
    message(crayon::green(paste0("Cross-validating supplied lme4 model:\n\t", 
                                 format(stats::formula(mod_obj$model)), "\n")))
  }
  if (inherits(mod_obj, "HibModelINLA")){
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
    message(crayon::green(paste0("Cross-validating supplied INLA model:\n\t", 
                                 m_form, "\n")))
    message(crayon::green("Have some patience, this may take a while...\n"))
  } 
  
  if (inherits(mod_obj, "HibModelBRMS")){
    message(crayon::green(paste0("Cross-validating supplied brms model:\n\t", 
                                 format(mod_obj$model$formula$formula), "\n")))
    
    message(crayon::green("Have some patience, this may take a while...\n"))
  }

  mod_specs <- get_model_info(mod_obj, re_form, use_original_mesh)
  # print(mod_specs)
  # Update within-function vars based on mod_specs
  rand_vars <- mod_specs$rand_vars
  re_form_vars <- mod_specs$re_form_vars
  re_form <- mod_specs$re_form

  # Random effects to consider - re_form - message?
  # paste(as.character(re_form), collapse = "; ")

  cv_out <- purrr::pmap(cv_spec,
                        function(fold_type, fold_col,  
                                 mod_ = mod_obj,
                                 db_ = db_obj,
                                 rand_seed = random_seed,
                                 re_f = re_form,
                                 backtr = backtransform,
                                 # orig_mesh = use_original_mesh,
                                 mod_specs_ = mod_specs){
                          hib_cv_mod_inner(mod_, db_, fold_type, fold_col, 
                                           rand_seed, re_f, backtr, #orig_mesh,
                                           mod_specs_)
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
  

  # Useful helpers - inla_fit_pred?

  insample_df_ls <- vector("list", length(re_form))

  if (inherits(mod_obj, "HibModelLME4")){
    # In-sample preds
    # orig_form <- stats::formula(mod_obj$model)
    orig_form <- mod_specs$lme4_specs$orig_form

    for (i in seq_along(re_form)){
      tmp_df <- mod_obj$data
      tmp_df$Random_effect_form <- as.character(re_form[i])
      insample_df_ls[[i]] <- tmp_df
    }
    
    insample_df <- predict_lme4_brms_cv(mod_obj$model, insample_df_ls, re_form)
  }
  
  
  if (inherits(mod_obj, "HibModelINLA")){
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # message(crayon::yellow("Warning: INLA insample prediction in development, current output is preliminary!"))
    
    # In-sample preds
    # insample_df <- mod_obj$data
    insample_df_ls <- vector("list", length(re_form))

    # Need orig_form (user supplied$response_fixed) to edit/backtransform response 
    orig_form <- mod_obj$user_formula$response_fixed
    fixed_covars <- all.vars(orig_form[3])
    fixed_covars <- fixed_covars[!fixed_covars == "Intercept"]
    
    fitted_form_ <- mod_obj$inla_fitted_formula
    orig_fam <- mod_obj$model$.args$family

    # ??: use mod_specs$inla_specs
    # orig_form <- mod_specs$inla_specs$orig_form
    # fixed_covars <- all.vars(orig_form[3])
    # fitted_form_ <- mod_specs$inla_specs$fitted_form_
    # orig_fam <- mod_specs$inla_specs$orig_fam
    
    stk_te_ls <- vector("list", length(re_form))
    # Make a new inla stack for test predictions
    stk_te <- mod_obj$inla_stack
    # Drop y, rename index
    stk_te$data$data$y <- NA
    
    # loop through re_form, set col to na if not in re_form_vars, but in rand_vars
    for (i in seq_along(re_form)){
      rand_vars_NA <- rand_vars[!rand_vars %in% re_form_vars[[i]]]
      tmp_stk <- stk_te
      names(tmp_stk$data$index)[names(tmp_stk$data$index) == "fit"] <- paste0("te__",
                                                                              paste(re_form_vars[[i]], collapse = "_"))
      if (length(rand_vars_NA) > 0){
        tmp_stk$effects$data <- tmp_stk$effects$data |> 
          dplyr::mutate_at(rand_vars_NA, ~ NA)

        stk_te_ls[[i]] <- tmp_stk
      } else {
        stk_te_ls[[i]] <- tmp_stk
      }

      # %%%%%%
      # TODO: is mod_obj$data same as tmp_stk$effects$data?
      # %%%%%%
      tmp_df <- mod_obj$data
      tmp_df$Random_effect_form <- as.character(re_form[i])
      insample_df_ls[[i]] <- tmp_df
    }

      # X_te <- dplyr::bind_rows(pred_df_ls)
    
    
    # stk_te_ls
    # # Drop iid random effects
    # # rand_vars <- all.vars(mod_obj$user_formula$random)
    
    # # TODO: what about incl some based on a re_form arg??
    # # for (i in seq_along(re_form){stk_te_i <- ...}

    # stk_te$effects$data <- stk_te$effects$data |> 
    #   dplyr::mutate_at(rand_vars, ~ NA)

    # print(head(stk_te$effects$data))
    # Combine stacks - does this work?
    comb_stk <- INLA::inla.stack(mod_obj$inla_stack,
                                 stk_te_ls[[1]])
    if (length(seq_along(re_form)) > 1){
      for (i in 2:length(re_form)){
        comb_stk <- INLA::inla.stack(comb_stk,
                                   stk_te_ls[[i]])
      }
    }
    
    # comb_stk <- INLA::inla.stack(mod_obj$inla_stack,
    #                              unlist(stk_te_ls))
    
    if (length(fixed_covars) > 0){
      check_flag <- 0
      while(check_flag == 0){
        tmp_m <- INLA::inla(formula = fitted_form_, # orig_form,
                            family = orig_fam,
                            data = INLA::inla.stack.data(comb_stk),
                            control.fixed = list(expand.factor.strategy='inla'),
                            # control.inla = list(strategy = "", int.strategy = ""), 
                            control.compute = list(dic = TRUE, waic = TRUE),
                            control.predictor = list(A = INLA::inla.stack.A(comb_stk),
                                                      link = 1),
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
                          control.predictor = list(A = INLA::inla.stack.A(comb_stk),
                                                    link = 1),
                          verbose = F)
    }
    
    
    # TODO
    for (i in seq_along(re_form)){
      tag_ <- paste0("te__", paste(re_form_vars[[i]], collapse = "_"))


      # Extract prediction and cis to insample_df
      pred_idx <- INLA::inla.stack.index(comb_stk, tag = tag_)$data
      tmp_pred <- tmp_m$summary.fitted.values[pred_idx,]
      
      insample_df_ls[[i]]$prediction <- tmp_pred$mean
      insample_df_ls[[i]]$prediction_q2.5 <- tmp_pred$`0.025quant`
      insample_df_ls[[i]]$prediction_q97.5 <- tmp_pred$`0.975quant`
      # Drop y - will be regenerated below for error estimation etc
      insample_df_ls[[i]] <- insample_df_ls[[i]] |> 
        dplyr::select(-y)
    }

    insample_df <- dplyr::bind_rows(insample_df_ls)
    # %%%%%
    # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  }
  
  if (inherits(mod_obj, "HibModelBRMS")){
    # In-sample preds
    orig_form <- mod_obj$model$formula$formula
 
    for (i in seq_along(re_form)){
      tmp_df <- mod_obj$data
      tmp_df$Random_effect_form <- as.character(re_form[i])
      insample_df_ls[[i]] <- tmp_df
    }
    
    insample_df <- predict_lme4_brms_cv(mod_obj$model, insample_df_ls, re_form)
    
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
  
  insample_metrics <- purrr::map(as.character(re_form),
            function(x, df = insample_df){
              df_sub <- df |>
                dplyr::filter(Random_effect_form == x)
              tmp <- calc_cv_metrics(df_sub$prediction, 
                                    df_sub$y, 
                                    df_sub$error)
              tmp$Random_effect_form <- x
              tmp
              }) |> 
    dplyr::bind_rows()
  # insample_metrics <- calc_cv_metrics(insample_df$prediction, 
  #                                     insample_df$y, 
  #                                     insample_df$error)
  
  insample_res = list(insample_pred_df = insample_df, 
                      insample_metrics = insample_metrics)
  
  # %%%%%
  
  # %%%
  # Sort out use_original_mesh
  # if not an inla model, set to NA
  if (!inherits(mod_obj, "HibModelINLA")){
    use_original_mesh <- NA
  }

  # if an inla, and spatial model, check it's not null, if it is, NA, else logical T/F
  if (inherits(mod_obj, "HibModelINLA")){
    if (mod_obj$inla_incl_spatial_term == TRUE){
      use_original_mesh <- mod_specs$inla_specs$use_original_mesh
    }
    # if an inla, but non spatial model, set to na
    if (mod_obj$inla_incl_spatial_term == FALSE){
      use_original_mesh <- NA
    }
  }
  # %%%
  
  # Return a cv object.. - incl original model? + df of cv predictions 
  # initialise cv_output...
  message(crayon::green("\nCollating cross-validation output ..."))
  
  cv_obj <- HibCV$new(
    mod_obj = mod_obj,
    cv_res = cv_out,
    insamp_res = insample_res,
    re_form = re_form,
    use_orig_mesh = use_original_mesh)
  
  return(cv_obj)
}


#' Cross-validate a fitted hib_XX_model
#' 
#' @description 
#' This function re-fits the model on different subsets of the data, so that subsequent 
#' evaluation of model performance can take place.
#' @details 
#' Expects a HibModel object.
#' @param mod_obj One of; [`HibModelLME4`], [`HibModelINLA`] or [`HibModelBRMS`].
#' @param db_obj A [`HibDatabase`] object, or, if modelling used multiple data-
#' sources, a list of database objects
#' @param fold_type A [`character`] indicating the type of cross-validation; one
#' of: "Random-2", "Random-3", "Random-5", "Random-10", "Source_ID", "SS", "SSBS",
#' "Country", "UN_subregion", "UN_region", "Biome", "Realm", "Custom".
#' @param fold_col A [`character`] indicating the column in the model data to 
#' use to define cross-validation folds
#' @param random_seed A [`numeric`] indicating the random seed used to define
#' random cross-validation folds - only used if 'Random-X' in fold_type
#' @param backtransform A function to convert the model prediction back to the original
#' scale of interest
#' @param re_form A [`list`] indicating the random effects variables to be
#' included in the cross-validation/insample predicitons - if multiple items supplied,
#' each option will be used in a separate cv prediction.
#' Values can be either NA (no random effects), NULL (all random effects), or a
#' rhs formula specifying particular random effects to include.
#' Note: only test set data points with random-effect values present in the 
#' training data will have predictions if re_form is not NA.
#' @param use_original_mesh A logical value indicating whether cross-validation
#' of spatial INLA models should use the original mesh for all cv folds. Computationally
#' this is slightly more efficient, but may affect cv performance metrics - witheld
#' sites will still exist in the mesh. - NO LONGER USED
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
                             re_form = list(NA),
                             backtransform = NULL,
                             # use_original_mesh = NULL, # TODO: drop - now in mod_specs
                             mod_specs = NULL){
  
  # Assertion checks 
  assertthat::assert_that(
    inherits(mod_obj, "HibModelLME4") || inherits(mod_obj, "HibModelINLA") || inherits(mod_obj, "HibModelBRMS"),
    inherits(db_obj, "HibDatabase") || ( is.list(db_obj) & all(unlist(lapply(db_obj, function(x){inherits(x, "HibDatabase")}))) ),
    is.character(fold_type),
    is.character(fold_col),
    is.numeric(random_seed) || is.null(random_seed),
    is.list(re_form),
    is.function(backtransform) || is.null(backtransform),
    # is.logical(use_original_mesh) || is.null(use_original_mesh),
    is.list(mod_specs)
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
    if(inherits(db_obj, "HibDatabase")){
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
  
  
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Get important vars from mod_spec...
  orig_data <- mod_specs$orig_data
  rand_vars <- mod_specs$rand_vars
  re_form_vars <- mod_specs$re_form_vars
  fix_fact_vars <- mod_specs$fix_fact_vars

  # %%%%%%%%%%%%%%%%%%
  # Below is specific to each cv type... - i.e. needed per cv_inner
  
  
  # If model is based on contrasts (i.e., comp-sim style), need to get, info on both
  # SSBS1 and SSBS2... - below, am just using SSBS1...
  
  # Add important fold variable - from db_obj - based on SSBS, or, 
  # if SSBS1/SSBS2 in orig_data, use these...
  if ("SSBS" %in% colnames(orig_data)){
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
    
    return(list(cv_pred_df  = NULL, 
                cv_metrics  = NULL,
                cv_fold_metrics  = NULL,
                cv_coef_df  = NULL,
                fold_type   = fold_type,
                fold_col    = fold_col,
                n_folds     = n_folds,
                random_seed = NULL
    ))
  }
  
  message(crayon::green(paste0("\tTotal number of folds: ", n_folds, "\n\n\tProgress")))
  
  # TODO: this is only really good for LME4 models, unless we make brms/inla
  # silent
  if ( (inherits(mod_obj, "HibModelLME4")) | 
        (inherits(mod_obj, "HibModelINLA")) ){
    # print("WOOP - lme4/inla")
    pb <- utils::txtProgressBar(min = 1, max = max(2, n_folds), style = 3)
  }
  # Get orig_form as a local var..

  if (inherits(mod_obj, "HibModelLME4")){
    incl_weights <- mod_specs$lme4_specs$incl_weights
    orig_form <- mod_specs$lme4_specs$orig_form
  }
  if (inherits(mod_obj, "HibModelINLA")){
    # incl_weights <- mod_specs$lme4_specs$incl_weights
    incl_weights <- FALSE
    orig_form <- mod_specs$inla_specs$orig_form
    use_original_mesh <- mod_specs$inla_specs$use_original_mesh
  }
  if (inherits(mod_obj, "HibModelBRMS")){
    incl_weights <- mod_specs$brms_specs$incl_weights
    orig_form <- mod_specs$brms_specs$orig_form
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


    if (incl_weights == TRUE){
      if (inherits(mod_obj, "HibModelLME4")){
        # TODO: this weighting is approapriate for lme4, but I don't think it is
        # correct for brms - can leave brms as is?
        weight_df <- orig_data |>
          dplyr::group_by(dataset_modelling_id) |>
          dplyr::summarise(prop_weight = sum(data_weight))
        
        tr <- tr |>
          dplyr::select(-data_weight)
        
        tr <- dplyr::left_join(tr, weight_df) |>
          dplyr::group_by(dataset_modelling_id) |>
          dplyr::mutate(data_weight = prop_weight/dplyr::n())
      }
      te <- te |>
          dplyr::select(-data_weight)
      
    }
    
    # Given each re_form item, make an output df for each
    pred_df_ls <- vector("list", length(re_form))
    for (j in seq_along(re_form)){
      # Filter 
      if (is.na(re_form_vars[[j]])){
        tmp <- te
      }
      if (is.character(re_form_vars[[j]])){
        tmp <- te
        for (var in re_form_vars[[j]]){
          tmp <- tmp |> 
            dplyr::filter(!!dplyr::sym(var) %in% tr[[var]])
        }
      }

      # TODO: need to check that factor vars in te are present in tr...
      for (var in fix_fact_vars){
        # print(var)
        tmp <- tmp |> 
            dplyr::filter(!!dplyr::sym(var) %in% tr[[var]])
      }
      pred_df_ls[[j]] <- tmp |> droplevels()
      # TODO: need to check that pred_df_ls[[j]] has nrow() > 0
      # If not... it can't be used...
      if (nrow(pred_df_ls[[j]]) > 0){
        pred_df_ls[[j]]$Random_effect_form <- as.character(re_form[j])  
      }
      # print(nrow(pred_df_ls[[i]]))
      # print("AD1_2004__Quaranta 1" %in% tr$SS)
    }
    # Note: sometimes, cv may not be feasible - e.g. SS cv, with an SS re_form - no data in te...



    if (inherits(mod_obj, "HibModelLME4")){
      # %%%%%%%%%%
      # LME4-specific
      # %%%%%%%%%%
      # fit model 
      # - lmer v glmer - # - family
      if (is.null(mod_specs$lme4_specs$fam)){
        if (mod_specs$lme4_specs$incl_weights == FALSE){
          tmp_m <- lme4::lmer(formula = mod_specs$lme4_specs$orig_form, 
                              data = tr,
                              REML = mod_specs$lme4_specs$REML)
        }
        if (mod_specs$lme4_specs$incl_weights == TRUE){
          tmp_m <- lme4::lmer(formula = mod_specs$lme4_specs$orig_form, 
                              data = tr,
                              REML = mod_specs$lme4_specs$REML,
                              weights = data_weight)
        }
      }
      if (!is.null(mod_specs$lme4_specs$fam)){
        if (mod_specs$lme4_specs$incl_weights == FALSE){
          tmp_m <- lme4::glmer(formula = mod_specs$lme4_specs$orig_form, 
                               data = tr,
                               family = mod_specs$lme4_specs$fam)
        }
        if (mod_specs$lme4_specs$incl_weights == TRUE){
          tmp_m <- lme4::glmer(formula = mod_specs$lme4_specs$orig_form, 
                               data = tr,
                               family = mod_specs$lme4_specs$fam,
                               weights = data_weight)
        }
      }
      
      te <- predict_lme4_brms_cv(tmp_m, pred_df_ls, re_form)
      
      # get coefs, and cis...
      fixef_df <- as.data.frame(lme4::fixef(tmp_m))
      colnames(fixef_df) <- "Coefficient_value"
      fixef_df <- tibble::rownames_to_column(fixef_df, "Coefficient_name")
      
      if (!is.null(mod_specs$lme4_specs$ci_type)){
        if (!is.na(mod_specs$lme4_specs$ci_type)){
          confint_df <- lme4::confint.merMod(tmp_m, parm = "beta_", method = mod_specs$lme4_specs$ci_type)
          confint_df <- as.data.frame(confint_df)
          colnames(confint_df) <- c("loCI", "hiCI")
          confint_df$Coefficient_name <- rownames(confint_df)
          fixef_df <- dplyr::left_join(fixef_df,
                                       confint_df,
                                       by = "Coefficient_name")  
        }
      }
      # %%%%%%%%%%
      # LME4-specific - END
      # %%%%%%%%%%
    }
    
    if (inherits(mod_obj, "HibModelINLA")){
      # %%%%%%%%%%
      # INLA-specific
      # %%%%%%%%%%
      
      # Make X_tr and X_te
      # print(head(tr))
      X_tr <- tr |> 
        dplyr::select(-y)
      
      # Drop iid random effects from te - where appropriate
      # rand_vars <- all.vars(orig_form_user_ls$random)
    
      # loop through pred_df_ls, set col to na if not in re_form_vars, but in rand_vars
      for (j in seq_along(re_form)){
        rand_vars_NA <- rand_vars[!rand_vars %in% re_form_vars[[j]]]

        if (length(rand_vars_NA) > 0){
          pred_df_ls[[j]] <- pred_df_ls[[j]] |> 
            dplyr::select(-y) |>
            dplyr::mutate_at(rand_vars_NA, ~ NA)
        } else {
          pred_df_ls[[j]] <- pred_df_ls[[j]] |> 
            dplyr::select(-y)
        }
        
      }

      X_te <- dplyr::bind_rows(pred_df_ls)
      
      # TODO: how to deal with no row X_te...

      # print(dim(X_tr))
      # print(dim(X_te))
      
      # If not a spatial model
      if (mod_specs$inla_specs$orig_incl_spatial_term == FALSE){
        # Just use model dfs
        stk_tr <- INLA::inla.stack(tag = "tr",
                                   data = list(y = tr$y),
                                   A = list(1),
                                   effects = list(X = X_tr))
        
        if (nrow(X_te) == 0){
          stk_te <- NULL  
        } else {
          stk_te <- INLA::inla.stack(tag = "te",
                                     data = list(y = NA),
                                     A = list(1),
                                     effects = list(X = X_te))  
        }
        
      }
      
      # TODO:
      # If spatial, make and generate mesh...
      if (mod_specs$inla_specs$orig_incl_spatial_term == TRUE){
        
        # print(orig_mesh_settings)
        # print()
        # Get lat/long
        if (("Longitude" %in% colnames(X_tr)) &
            ("Latitude" %in% colnames(X_tr)) ){
          coords_tr <- as.matrix(X_tr[,c("Longitude", "Latitude")])
          if (nrow(X_te) == 0){
            coords_te <- NULL  
          } else {
            coords_te <- as.matrix(X_te[,c("Longitude", "Latitude")])
          }
        }
        if (("Longitude_midpoint" %in% colnames(X_tr)) &
            ("Latitude_midpoint" %in% colnames(X_tr)) ){
          coords_tr <- as.matrix(X_tr[,c("Longitude_midpoint", "Latitude_midpoint")])
          if (nrow(X_te) == 0){
            coords_te <- NULL  
          } else {
            coords_te <- as.matrix(X_te[,c("Longitude_midpoint", "Latitude_midpoint")])
          }
        }
        
        if (use_original_mesh == TRUE){
          cv_mesh <- mod_specs$inla_specs$orig_mesh
        }
        
        if (use_original_mesh == FALSE){
          if (mod_specs$inla_specs$mesh_source == "mesh_settings"){
            message(crayon::green("Building mesh from original settings.\n"))
            # Need to re-build a mesh - based only on tr points
            # Spatial settings
            # orig_mesh
            # orig_mesh_settings
            if (mod_specs$inla_specs$orig_mesh_settings$mesh_dim == "3d"){
              # print("Building 3d mesh...")
              coords_tr <- INLA::inla.mesh.map(coords_tr, projection = "longlat")
              cv_mesh <- INLA::inla.mesh.create(loc = coords_tr,
                                                   globe = 8)
            }
            
            if (mod_specs$inla_specs$orig_mesh_settings$mesh_dim == "2d"){
              # print("Building 2d mesh...")
              # TODO: Why is building the cv mesh taking soooo long...
              cv_mesh <- INLA::inla.mesh.2d(loc = coords_tr,
                                               max.edge = mod_specs$inla_specs$orig_mesh_settings$max.edge,
                                               cutoff = mod_specs$inla_specs$orig_mesh_settings$cutoff,
                                               offset = mod_specs$inla_specs$orig_mesh_settings$offset)
              # print("Mesh built!")
            }
          }
          
          if (mod_specs$inla_specs$mesh_source == "user_provided"){
            cv_mesh <- mod_specs$inla_specs$orig_mesh
            # TODO: CHECK THIS IS WORKING
            message(crayon::yellow("Warning: mesh in original model was user-provided, using this for cross-validation.\n"))
          }
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
        
        
        # Combine into stacks - tr and te
        stk_tr <- INLA::inla.stack(tag = "tr",
                                   data = list(y = tr$y),
                                   A = list(1, A_tr),
                                   effects = list(X = X_tr,
                                                  s = s.site))
        
        # -for te
        if (nrow(X_te) == 0){
          stk_te <- NULL  
        } else {
          A_te <- INLA::inla.spde.make.A(cv_mesh, 
                                         loc = coords_te)
          stk_te <- INLA::inla.stack(tag = "te",
                                     data = list(y = NA),
                                     A = list(1, A_te),
                                     effects = list(X = X_te,
                                                    s = s.site)) 
        }
      }
      
      # Combine stacks
      if (is.null(stk_te)){
        comb_stk <- stk_tr
      } else {
        comb_stk <- INLA::inla.stack(stk_tr, 
                                     stk_te)  
      }
      
      # Fit model
      fixed_covars <- all.vars(mod_specs$inla_specs$orig_form_user_ls$response_fixed[[3]])
      
      if (length(fixed_covars) > 0){
        check_flag <- 0
        while(check_flag == 0){
          tmp_m <- INLA::inla(formula = mod_specs$inla_specs$fitted_form_, # orig_form,
                           family = mod_specs$inla_specs$orig_fam,
                           data = INLA::inla.stack.data(comb_stk),
                           control.fixed = list(expand.factor.strategy='inla'),
                           # control.inla = list(strategy = "", int.strategy = ""), 
                           control.compute = list(dic = TRUE, waic = TRUE),
                           control.predictor = list(A = INLA::inla.stack.A(comb_stk)),
                           verbose = T)
          
          check_flag <- nrow(tmp_m$summary.fixed)
        }
      } else {
        tmp_m <- INLA::inla(formula = mod_specs$inla_specs$fitted_form_, # orig_form,
                         family = mod_specs$inla_specs$orig_fam,
                         data = INLA::inla.stack.data(comb_stk),
                         control.fixed = list(expand.factor.strategy='inla'),
                         # control.inla = list(strategy = "", int.strategy = ""), 
                         control.compute = list(dic = TRUE, waic = TRUE),
                         control.predictor = list(A = INLA::inla.stack.A(comb_stk)),
                         verbose = T)
      }
      
      
      # Extract prediction and cis to te (a df of data + pred)
      if (nrow(X_te) == 0){
        te <- X_te
      } else {
        pred_idx <- INLA::inla.stack.index(comb_stk, tag = "te")$data
        tmp_pred <- tmp_m$summary.fitted.values[pred_idx,]
        
        # X_te <- dplyr::bind_cols(tmp_pred[,c("mean", "0.025quant", "0.975quant")])
        X_te$prediction <- tmp_pred$mean
        X_te$prediction_q2.5 <- tmp_pred$`0.025quant`
        X_te$prediction_q97.5 <- tmp_pred$`0.975quant`
        
        if ("y" %in% colnames(X_te)){
          te <- X_te |> 
            dplyr::select(-y)
        } else {
          te <- X_te
        }
      }
      
      
      # # Drop y from te - will be regenerated below for error estimation etc
      # te <- te |> 
      #   dplyr::select(-y)
      
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
    
    if (inherits(mod_obj, "HibModelBRMS")){
      # %%%%%%%%%%
      # BRMS-specific
      # %%%%%%%%%%
      
      # fit model 
      tmp_args <- NULL
      tmp_args <- mod_specs$brms_specs$m_args
      tmp_args[["data"]] <- tr
      
      tmp_m <- do.call(brms::brm,
                       tmp_args)
      
      # predict or fitted??

      te <- predict_lme4_brms_cv(tmp_m, pred_df_ls, re_form)

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
    
    if ( (inherits(mod_obj, "HibModelLME4")) |
         (inherits(mod_obj, "HibModelINLA")) ){
      # print("PROGRESS UPDATE")
      utils::setTxtProgressBar(pb, i)
    }
  }
  
  if ( (inherits(mod_obj, "HibModelLME4")) |
        (inherits(mod_obj, "HibModelINLA")) ){
    close(pb)
  }
  
  cv_out_df <- dplyr::bind_rows(cv_out_ls)
  coef_out_df <- dplyr::bind_rows(coef_out_ls)
  
  # Collate cv outputs
  
  
  # Only return random seed if random-
  if (fold_type %in% c(c("Random-2", "Random-3", "Random-5", "Random-10" ))){
    rand_num <- random_seed
  } else {
    rand_num <- NA
  }
  
  # TODO: check nrow of cv_out_df, if 0, skip, return NULL for 
  if (nrow(cv_out_df) == 0){

    message(crayon::yellow(paste0("Warning: Due to combination of cross-validation approach ( ", 
                                  fold_type_msg, 
                                  " ) and re_form options specified ( ",
                                  as.character(re_form), 
                                  " ), no out-of-sample predictions have been made.\n")))
    
    return(list(cv_pred_df  = NULL, 
                cv_metrics  = NULL,
                cv_fold_metrics  = NULL,
                cv_coef_df  = coef_out_df,
                fold_type   = fold_type,
                fold_col    = fold_col,
                n_folds     = n_folds,
                random_seed = rand_num
    ))
  }
  
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
  # TODO:
  # per Random_effect_form/re_form
  cv_metrics <- purrr::map(as.character(re_form),
            function(x, df = cv_out_df, fld_typ = fold_type_msg){
              df_sub <- df |>
                dplyr::filter(Random_effect_form == x)
              # print(head(df_sub))
              if (nrow(df_sub) == 0){
                # TODO: message
                
                message(crayon::yellow(paste0("Warning: No predictions made for combination of '",
                                              fld_typ, "' cross-validation and random-effects structure '",
                                              x, "'.\n"
                                              )))
                return(NULL)
              }
              tmp <- calc_cv_metrics(df_sub$prediction, 
                                    df_sub$y, 
                                    df_sub$error)
              tmp$Random_effect_form <- x
              tmp
            }) |> 
    dplyr::bind_rows()
  cv_fold_metrics <- purrr::map(as.character(re_form),
            function(x, df = cv_out_df){
              df_sub <- df |>
                dplyr::filter(Random_effect_form == x)
              if (nrow(df_sub) == 0){
                return(NULL)
              }
              tmp <- calc_per_fold_cv_metrics(df_sub$prediction, 
                                    df_sub$y, 
                                    df_sub$error,
                                    df_sub$fold_id)
              tmp$Random_effect_form <- x
              tmp
            }) |> 
    dplyr::bind_rows()

  # cv_metrics <- calc_cv_metrics(cv_out_df$prediction, 
  #                               cv_out_df$y, 
  #                               cv_out_df$error)
  
  # cv_fold_metrics <- calc_per_fold_cv_metrics(cv_out_df$prediction, 
  #                                             cv_out_df$y, 
  #                                             cv_out_df$error,
  #                                             cv_out_df$fold_id)
  
  
  
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

  # NOTE: relative metrics here use mean of obs or similar as a benchmark
  # but this is the mean over all data and all cv-folds... - may be more appropriate to calc 
  # relative to 'in-sample' mean for data used to fit model...
  # Alternatively, if we just think of the denominator as a scaling factor to
  # align errors from various response scales, this is less of a concern - 
  # interpretation is just error scaled by average obs_minus_obs_mean or similar 
  # < 1 if better than average as a predictor, > 1 if worse.

  pred_mean <- mean(pred, na.rm = T)
  obs_mean <- mean(obs, na.rm = T)
  
  obs_minus_obs_mean <- obs - obs_mean
  
  # cor, corsq, 
  # Need to rm NAs...
  cor_out <- stats::cor(pred, obs, use = "complete.obs")
  cor_sq <- cor_out^2
  
  # me, mae, rae, mape, mse, rmse, nrmse, rmsle, rse, rrmse
  me <- mean(error, na.rm = T)
  mede <- stats::median(error, na.rm = T)

  # TODO:
  # A scaled bias measure - e.g. me[di]an( error / obs_minus_obs_mean, na.rm = T )
  # or m[ed]e / me[di]an(obs_minus_obs_mean, na.rm = T)

  ae <- abs(error)
  
  mae <- mean(ae, na.rm = T)
  
  # %%%%%
  # 20240214: Update included metrics, reduce redunctancy and add some new ones.
  # medae
  medae <- stats::median(ae, na.rm = T)

  # nmae
  # # nmae <- mae / obs_mean

  # rae (equivalent to mase - as formulted here) - done
  rae <- sum(ae, na.rm = T) / sum(abs(obs_minus_obs_mean), na.rm = T)

  # mrae (mean)
  mrae <- mean( ae / abs(obs_minus_obs_mean) , na.rm = T )

  # medrae (median)
  medrae <- stats::median( ae / abs(obs_minus_obs_mean), na.rm = T)

  # mape - done
  # smape - done

  # mse - done
  mse <- mean(error^2, na.rm = T)
  
  # nmse 
  nmse <- mse / stats::var(obs, na.rm = T)

  # rmse - done
  rmse <- sqrt(mse)
  
  # nrmse (mean) - done
  nrmse <- rmse / obs_mean
  # nrmse (sd)
  nrmse_sd <- rmse / stats::sd(obs, na.rm = T)

  # rrmse (rmse/rmse_naive)
  rrmse <- rmse / sqrt(mean(obs_minus_obs_mean^2 , na.rm = T))

  # logrrmse 
  logrrmse <- log(rrmse)
  
  # rmsle - done
  # rse - done
  # r2 - done

  # gmrae - done
  # %%%%%

  # are <- ??
  
  # rmae <- rae/n - or same? - should be the same, both divided by N.
  # mean of rae?: mrae
  # but, mean relative absolute error is slightly different...
  # mrae <- mean( ae / abs(obs_minus_obs_mean) , na.rm = T )

  # mape - issue of zeros in obs... - v small adjustment
  # issues of negative obs? - 
  # mape <- mean( (ae / obs), na.rm = T) * 100 OLD version
  # mape <- mean( (ae / abs(obs) ), na.rm = T) * 100 OLD version
  mape <- mean( abs( error / obs), na.rm = T ) * 100

  # smape <- mean( (ae / ( (pred + obs) / 2) ), na.rm = T ) * 100 OLD version
  # TODO: update smape to 0-100 scale
  smape <- mean( (ae / ( abs(pred) + abs(obs) ) ), na.rm = T ) * 100
  
  
  # nrmse_sd <- rmse / sd(pred, na.rm = T) OLD version - but maybe better?!
  # Alternative
  # nrmse <- rmse / obs_mean
  # ?? mean(abs(obs), na.rm = T) ??

  # rmse relative to rmse of mean - use this
  # rrmse <- rmse / sqrt(mean(obs_minus_obs_mean^2 , na.rm = T)) - same as sqrt(rse)??

  # # rrmse <- sqrt( ( sum(error^2, na.rm = T) / sum(pred^2, na.rm = T) ) ) OLD version
  # Check this
  # Alternative
  # rrmse <- rmse / obs_mean * 100

  rmsle <- sqrt(mean( ( log(obs + 1) - log(pred + 1) )^2 , na.rm = T))
  
  rse <- sum(error^2, na.rm = T) / sum( obs_minus_obs_mean^2 , na.rm = T)
  # This is a bit different to sum(error^2 / obs_minus_obs_mean^2) - but still okay??

  # Can sqrt rse to get back on scale of obs
  # - rrse 

  # rsq - [1 - (PRESS / sums of squares total)] * 100
  rsq <- ( 1 - rse )
  
  # TODO:
  # Any metric that is relative to "obs_mean" or similar could also be calculated per fold 
  # - global mean is not per-fold training sample mean...
  # rsq, rse, rae - all relative to using mean as a predictor...
  # need fold_id...
  
  # TODO:
  # median relative absolute error
  # mrae <- stats::median( abs( error / obs_minus_obs_mean ), na.rm = T)
  # median and mean versions?

  # geometric mean relative absolute error
  gmrae <- exp(mean(log( abs( error / obs_minus_obs_mean ) ), na.rm = T))
  
  # mean absolute scaled error: 
  # mase = mae / mean(abs(obs_minus_obs_mean), na.rm = T) # clearer.
  # mase <- mean( ae / mean(abs(obs_minus_obs_mean), na.rm = T), na.rm = T) # - equivalent
  # isn't this just rae - or at least v similar
  # similar, mase is sometime refered to as relative mean absolute error (mae/mae_naive)
  # I think this is just rae - but scaled by N, which cancels out...

  # Notmalised metrics:
  # Normalised MAE - denominator is obs_mean
  # nmae <- mae / obs_mean [ or mean of abs values: mean(abs(obs), na.rm = T)]
  # Normalised MedAE - ?? denominator as median of obs ??
  # Normalised RMSLE


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
  
  N_p <- length(pred)
  N_o <- length(pred)

  N_p_no_NA <- sum(!is.na(pred))

  # Make metr into a df
  
  # TO ADD
  # "Median_absolute_error"
  # "Mean_relative_absolute_error"
  # "Normalised_mean_squared_error"
  # "Normalised_root_mean_squared_error_sd"
  # "Relative_root_mean_squared_error"
  # "Log_relative_root_mean_squared_error"

  # TO DROP
  # "Mean_absolute_scaled_error"

  metr_df <- data.frame(Metric = c("Correlation",
                                   "Correlation_squared",
                                   "Mean_error", 
                                   "Median_error", 
                                   "Mean_absolute_error", 
                                   "Median_absolute_error", # New
                                   "Relative_absolute_error", 
                                   "Mean_relative_absolute_error", # New
                                   "Median_relative_absolute_error", # Reposition
                                   "Mean_absolute_percentage_error", 
                                   "Symmetric_mean_absolute_percentage_error", 
                                   "Mean_squared_error", 
                                   "Normalised_mean_squared_error", # New
                                   "Root_mean_squared_error", 
                                   "Normalised_root_mean_squared_error",
                                   "Normalised_root_mean_squared_error_sd", # New
                                   "Relative_root_mean_squared_error", # Recalc
                                   "Log_relative_root_mean_squared_error", # New
                                   "Root_mean_squared_log_error",
                                   "Relative_squared_error", 
                                   
                                   "Geometric_mean_relative_absolute_error",
                                   # "Mean_absolute_scaled_error", # Drop
                                   "Predictive_R_squared"),
                        Metric_val = c(cor_out, cor_sq,
                                       me, mede, mae, medae, rae, mrae, medrae,
                                       mape, smape, 
                                       mse, nmse, rmse, nrmse, nrmse_sd, 
                                       rrmse, logrrmse, rmsle, rse, 
                                       gmrae,
                                       # mase,
                                       rsq),
                        Metric_lo = c(NA, NA, 
                                      me_lo, NA, NA, NA, NA, NA, NA, NA, NA,
                                      mse_lo, NA, rmse_lo, NA, NA,
                                      NA, NA, NA, NA, 
                                      NA, 
                                      # NA, 
                                      NA),
                        Metric_hi = c(NA, NA, 
                                      me_hi, NA, NA, NA, NA, NA, NA, NA, NA,
                                      mse_hi, NA, rmse_hi, NA, NA,
                                      NA, NA, NA, NA, 
                                      NA, 
                                      # NA, 
                                      NA),
                        Metric_shrt = c("cor", "corsq",
                                        "me", "mede", "mae", "medae", "rae", "mrae", "medrae",
                                        "mape", "smape", 
                                        "mse", "nmse", "rmse", "nrmse", "nrmse_sd",
                                        "rrmse", "lrrmse", "rmsle", "rse", 
                                        # "mrae",
                                        "gmrae",
                                        # "mase",
                                        "rsq")) |>
    dplyr::mutate(Metric_fct = factor(Metric,
                                      levels = c("Correlation",
                                   "Correlation_squared",
                                   "Mean_error", 
                                   "Median_error", 
                                   "Mean_absolute_error", 
                                   "Median_absolute_error", # New
                                   "Relative_absolute_error", 
                                   "Mean_relative_absolute_error", # New
                                   "Median_relative_absolute_error", # Reposition
                                   "Mean_absolute_percentage_error", 
                                   "Symmetric_mean_absolute_percentage_error", 
                                   "Mean_squared_error", 
                                   "Normalised_mean_squared_error", # New
                                   "Root_mean_squared_error", 
                                   "Normalised_root_mean_squared_error",
                                   "Normalised_root_mean_squared_error_sd", # New
                                   "Relative_root_mean_squared_error", # Recalc
                                   "Log_relative_root_mean_squared_error", # New
                                   "Root_mean_squared_log_error",
                                   "Relative_squared_error", 
                                   
                                   "Geometric_mean_relative_absolute_error",
                                   # "Mean_absolute_scaled_error", # Drop
                                   "Predictive_R_squared"),
                                      labels = c("Correlation",
                                                 "Correlation^2",
                                                 "ME~(Bias)",
                                                 "Median~error", 
                                                 "MAE", 
                                                 "MedAE",
                                                 "RAE", 
                                                 "MRAE",
                                                 "MedRAE",
                                                 "MAPE", 
                                                 "sMAPE", 
                                                 "MSE", 
                                                 "NMSE",
                                                 "RMSE", 
                                                 "NRMSE",
                                                 "NRMSE~(SD)",
                                                 "RRMSE",
                                                 "LRRMSE",
                                                 "RMSLE",
                                                 "RSE", 
                                                 # "MRAE",
                                                 "GMRAE",
                                                 # "MASE",
                                                 "Predictive~R^2")),
      N = N_p_no_NA)
  
  
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

                     mede = stats::median(error, na.rm = T),
                     mae = mean(abs(error), na.rm = T),
                     medae = stats::median(abs(error), na.rm = T),
                     
                     # TOCHECK
                     rae = sum(abs(error), na.rm = T) / sum(abs(obs - obs_mean), na.rm = T),
                     mrae = mean( abs(error) / abs(obs - obs_mean), na.rm = T),
                     medrae = stats::median( abs(error) / abs(obs - obs_mean), na.rm = T),

                     # mape - issue of zeros in obs... - v small adjustment
                     mape = mean( (abs(error) / obs), na.rm = T) * 100, 
                     # smape = mean( (abs(error) / ( (pred + obs) / 2) ), na.rm = T ) * 100,OLD version
                     smape = mean( (abs(error) / ( abs(pred) + abs(obs) ) ), na.rm = T ) * 100,
                     
                     mse = mean(error^2, na.rm = T),
                     nmse = mean(error^2, na.rm = T) / stats::var(obs, na.rm = T),
                     
                     rmse = sqrt(mse),
                     # nrmse = rmse / sd(pred, na.rm = T),OLD version
                     nrmse = rmse / obs_mean,
                     # Multiple nrmse options - scaled by sd, mean, range, quartile range
                     nrmse_sd = rmse / stats::sd(obs, na.rm = T),

                     rrmse = rmse / sqrt(mean( (obs - obs_mean)^2 , na.rm = T)),
                     logrrmse = log(rrmse),
                     # rrmse = sqrt( ( sum(error^2, na.rm = T) / sum(pred^2, na.rm = T) ) ), OLD version
                     # rrmse = rmse / obs_mean * 100, # OLD
                     rmsle = sqrt(mean( ( log(obs + 1) - log(pred + 1) )^2 , na.rm = T)),
                     # %%%
                     
                     rse = sum(error^2, na.rm = T) / sum( (obs - obs_mean)^2 , na.rm = T),

                     # geometric mean relative absolute error
                     gmrae = exp(mean(log( abs( error / (obs - obs_mean) ) ), na.rm = T)),
                    
                     # mean absolute scaled error: 
                     # mase = mean( abs(error) / mean(abs((obs - obs_mean)), na.rm = T), na.rm = T), #- equivalent

                     rsq = 1 - rse ,

                     N = sum(!is.na(pred))
                     ) |>
    dplyr::select(-c(pred_mean, obs_mean)) |>
    # dplyr::select(fold_id, N, cor_out, cor_sq, me, mede, mae, mape, smape, mse, rmse, nrmse,
    #               rrmse, rmsle,
    #               rae, rse, rsq) |> # N
    tidyr::pivot_longer(!c(fold_id, N), 
                        names_to = "Metric_shrt", 
                        values_to = "Metric_val") |>
    dplyr::mutate(Metric_shrt = dplyr::case_when(Metric_shrt == "cor_out" ~ "cor",
                                                  Metric_shrt == "cor_sq" ~ "corsq",
                                                  T ~ Metric_shrt)) |>
    dplyr::mutate(Metric = dplyr::case_when(Metric_shrt == "cor"   ~ "Correlation",
                                            Metric_shrt == "corsq" ~ "Correlation_squared",
                                            Metric_shrt == "me"    ~ "Mean_error", 
                                            Metric_shrt == "mede"  ~ "Median_error", 
                                            Metric_shrt == "mae"   ~ "Mean_absolute_error", 
                                            Metric_shrt == "medae"  ~ "Median_absolute_error", # New
                                            Metric_shrt == "rae"   ~ "Relative_absolute_error", 
                                            Metric_shrt == "mrae"  ~ "Mean_relative_absolute_error",  # New
                                            Metric_shrt == "medrae"~ "Median_relative_absolute_error", # Reposition
                                            Metric_shrt == "mape"  ~ "Mean_absolute_percentage_error", 
                                            Metric_shrt == "smape" ~ "Symmetric_mean_absolute_percentage_error", 
                                            Metric_shrt == "mse"   ~ "Mean_squared_error", 
                                            Metric_shrt == "nmse"  ~ "Normalised_mean_squared_error", # New
                                            Metric_shrt == "rmse"  ~ "Root_mean_squared_error", 
                                            Metric_shrt == "nrmse" ~ "Normalised_root_mean_squared_error",
                                            Metric_shrt == "nrmse_sd"~ "Normalised_root_mean_squared_error_sd", # New
                                            Metric_shrt == "rrmse" ~ "Relative_root_mean_squared_error",
                                            Metric_shrt == "lrrmse"~ "Log_relative_root_mean_squared_error", # New
                                            Metric_shrt == "rmsle" ~ "Root_mean_squared_log_error",
                                            Metric_shrt == "rse"   ~ "Relative_squared_error", 

                                            # Metric_shrt == "mrae"  ~ "Median_relative_absolute_error", 
                                            Metric_shrt == "gmrae" ~ "Geometric_mean_relative_absolute_error", 
                                            # Metric_shrt == "mase"  ~ "Mean_absolute_scaled_error", Drop
                                            
                                            Metric_shrt == "rsq"   ~ "Predictive_R_squared",
                                            T ~ NA),
                  Metric_fct = factor(Metric,
                                      levels = c("Correlation",
                                   "Correlation_squared",
                                   "Mean_error", 
                                   "Median_error", 
                                   "Mean_absolute_error", 
                                   "Median_absolute_error", # New
                                   "Relative_absolute_error", 
                                   "Mean_relative_absolute_error", # New
                                   "Median_relative_absolute_error", # Reposition
                                   "Mean_absolute_percentage_error", 
                                   "Symmetric_mean_absolute_percentage_error", 
                                   "Mean_squared_error", 
                                   "Normalised_mean_squared_error", # New
                                   "Root_mean_squared_error", 
                                   "Normalised_root_mean_squared_error",
                                   "Normalised_root_mean_squared_error_sd", # New
                                   "Relative_root_mean_squared_error", # Recalc
                                   "Log_relative_root_mean_squared_error", # New
                                   "Root_mean_squared_log_error",
                                   "Relative_squared_error", 
                                   
                                   "Geometric_mean_relative_absolute_error",
                                   # "Mean_absolute_scaled_error", # Drop
                                   "Predictive_R_squared"),
                                      labels = c("Correlation",
                                                 "Correlation^2",
                                                 "ME~(Bias)",
                                                 "Median~error", 
                                                 "MAE", 
                                                 "MedAE",
                                                 "RAE", 
                                                 "MRAE",
                                                 "MedRAE",
                                                 "MAPE", 
                                                 "sMAPE", 
                                                 "MSE", 
                                                 "NMSE",
                                                 "RMSE", 
                                                 "NRMSE",
                                                 "NRMSE~(SD)",
                                                 "RRMSE",
                                                 "LRRMSE",
                                                 "RMSLE",
                                                 "RSE", 
                                                 # "MRAE",
                                                 "GMRAE",
                                                 # "MASE",
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




#' Helper function to collate information about fitted model
#' 
#' @description 
#' This function returns a list of model-specific information
#' @details 
#' Expects numeric inputs
#' @param mod_obj One of; [`HibModelLME4`], [`HibModelINLA`] or [`HibModelBRMS`].
#' @param re_form A [`list`] indicating the random effects variables to be
#' included in the cross-validation/insample predicitons.
#' @param use_original_mesh A logical value indicating whether cross-validation
#' of spatial INLA models should use the original mesh for all cv folds.
#' @keywords model, evaluation, processing
#' @examples 
#' \dontrun{
#'   model_info <- get_model_info(mod_lme4_obj, list(NA, ~ (1|SS)))
#' } 
#' @noRd
get_model_info <- function(mod_obj, re_form, use_original_mesh){
  
  # Assertion checks 
  assertthat::assert_that(
    inherits(mod_obj, "HibModelLME4") || inherits(mod_obj, "HibModelINLA") || inherits(mod_obj, "HibModelBRMS"),
    is.list(re_form),
    is.logical(use_original_mesh) || is.null(use_original_mesh)
  )
  
  # Should we also try to "keep" SS together in other blockings...
  if (inherits(mod_obj, "HibModelLME4")){
    rand_vars <- names(lme4::ranef(mod_obj$model))
  }
  if (inherits(mod_obj, "HibModelBRMS")){
    rand_vars <- names(brms::ranef(mod_obj$model))
  }
  if (inherits(mod_obj, "HibModelINLA")){
    rand_vars <- all.vars(mod_obj$user_formula$random)
  }
  
  re_form_vars <- purrr::map(re_form,
           function(x, rand_vars_ = rand_vars){
             if (is.null(x)){
               out <- rand_vars_
               return(out)
             }
             if (is.na(c(x))){
               out <- NA
             }
             if (inherits(x, "formula")){
               out <- all.vars(x)
               # Check all are in rand_vars 
               assertthat::assert_that(all(out %in% rand_vars),
                                      msg = "Specified random-effects structure contains variables not in the original model.")
             }
             out
           })
  
  # TODO: check for re_form_vars duplicates
  # How do drop... - just take the first...
  # Message user...
  re_form_usr <- re_form
  re_form <- re_form[!duplicated(re_form_vars)]
  re_form_vars <- re_form_vars[!duplicated(re_form_vars)]
  

  if (length(re_form) != length(re_form_usr)){
    message(crayon::yellow(paste0("Warning: duplication detected in specified ",
                                  "ranfom-effects structures.\n",
                                  "De-duplicated re_form: list(",
                                  paste(as.character(re_form), collapse = ", "), ")")))
  }


  mod <- mod_obj$model
  orig_data <- mod_obj$data
  
  # fixed_factor_vars
  # need to get... from fixeff
  

  incl_weights <- FALSE
  if ("data_weight" %in% colnames(orig_data)){
    incl_weights <- TRUE
  }
  
  lme4_specs <- NULL
  inla_specs <- NULL
  brms_specs <- NULL

  if (inherits(mod_obj, "HibModelLME4")){
    # %%%%%%%%%%
    # LME4-specific data prep 
    # %%%%%%%%%%
    orig_form <- stats::formula(mod)
    
    all_covars <- all.vars(stats::formula(stats::delete.response(stats::terms(orig_form))))
    # But what about random effects... - 
    rand_vars <- names(lme4::ranef(mod_obj$model))
    fix_vars <- all_covars[!all_covars %in% rand_vars]
    
    
    fam <- NULL
    if (inherits(mod@resp, "glmResp")){
      fam <- mod@resp$family
    }
    
    if (mod@devcomp$dims[["REML"]]){
      REML <- T
    } else {REML <- F}
    
    ci_type <- mod_obj$ci_type
    # %%%%%%%%%%
    # LME4-specific dataprep - END
    # %%%%%%%%%%

    lme4_specs <- list(orig_form    = orig_form,
                       fam          = fam,
                       REML         = REML,
                       ci_type      = ci_type,
                       incl_weights = incl_weights)
  }
  
  if (inherits(mod_obj, "HibModelINLA")){
    # %%%%%%%%%%
    # INLA-specific dataprep 
    # %%%%%%%%%%
    
    # Formulas 
    
    orig_form <- mod_obj$user_formula$full
    fitted_form <- stats::as.formula(mod$.args$formula)
    fitted_form_ <- mod_obj$inla_fitted_formula
    
    orig_form_user_ls <- mod_obj$user_formula
    # $full
    # $response_fixed
    # $random
    
    fix_vars <- all.vars(stats::formula(stats::delete.response(stats::terms(mod_obj$user_formula$response_fixed))))
    # ignore "Intercept" too
    fix_vars <- fix_vars[!fix_vars == "Intercept"]
    
    # family
    orig_fam <- mod$.args$family
    
    # Spatial settings
    orig_incl_spatial_term <- NA
    mesh_source <- NA
    orig_mesh <- NA
    orig_mesh_settings <- NA
    
    # print(use_original_mesh)
    orig_incl_spatial_term <- mod_obj$inla_incl_spatial_term
    if (orig_incl_spatial_term){
      if (is.null(use_original_mesh)){
        use_original_mesh <- TRUE
        message(crayon::yellow(paste0("Warning: use_original_mesh not specified, ",
                                      "using default value of TRUE.\n",
                                      "This should speed up computation, but may affect ",
                                      "estimates of model predictive performance.\n")))
      }
      # is mesh user provided or based on mesh_settings?
      mesh_source <- mod_obj$inla_mesh_source
      if (mesh_source == "user_provided"){
        orig_mesh <-  mod_obj$inla_mesh
      }
      if (mesh_source == "mesh_settings"){
        orig_mesh_settings <-  mod_obj$inla_mesh_settings
      }
      if (use_original_mesh == TRUE){
        orig_mesh <-  mod_obj$inla_mesh
        mesh_source <- "user_provided"
        message(crayon::green("Cross-validation will use spatial mesh from original model fit.\n"))
      }
      if (use_original_mesh == FALSE){
        message(crayon::green("Cross-validation will use a new spatial mesh for each fold.\n"))
      }
    }
    
    # prior_list ...
    # fit_settings ...
    # %%%%%%%%%%
    # INLA-specific dataprep -  END
    # %%%%%%%%%%
    inla_specs <- list(orig_form    = orig_form,
                       fitted_form  = fitted_form,
                       fitted_form_ = fitted_form_,
                       orig_form_user_ls = orig_form_user_ls,
                       orig_fam     = orig_fam,
                       orig_incl_spatial_term = orig_incl_spatial_term,
                       use_original_mesh = use_original_mesh,
                       mesh_source  = mesh_source, 
                       orig_mesh    = orig_mesh,
                       orig_mesh_settings = orig_mesh_settings)
  }
  
  if (inherits(mod_obj, "HibModelBRMS")){
    # %%%%%%%%%%
    # BRMS-specific dataprep -  
    # %%%%%%%%%%
    orig_form <- mod$formula$formula
    
    all_covars <- all.vars(stats::formula(stats::delete.response(stats::terms(orig_form))))
    # But what about random effects... - 
    rand_vars <- names(brms::ranef(mod_obj$model))
    fix_vars <- all_covars[!all_covars %in% rand_vars]
    
    orig_fam <- mod$family
    orig_extra_args <- mod_obj$extra_args
    
    m_args <- orig_extra_args
    m_args[["formula"]] <- orig_form
    m_args[["family"]] <- orig_fam
    # %%%%%%%%%%
    # BRMS-specific dataprep -  END
    # %%%%%%%%%%
    brms_specs = list(orig_form = orig_form,
                      orig_fam  = orig_fam,
                      orig_extra_args = orig_extra_args,
                      incl_weights = incl_weights,
                      m_args = m_args)
  }
  
  fix_fact_vars_idx <- purrr::map(fix_vars,
                              function(x, df = orig_data){
                                (is.factor(df[[x]]) | is.character(df[[x]]))
                              }) |>
    unlist()
  fix_fact_vars <- fix_vars[fix_fact_vars_idx]
  
  # print(fix_vars)
  # print(fix_fact_vars)
  
  return(list(mod = mod,
              orig_data = orig_data,
              
              fix_fact_vars = fix_fact_vars,
  
              rand_vars = rand_vars,
              re_form_vars = re_form_vars,
              re_form   = re_form,

              lme4_specs = lme4_specs,
              inla_specs = inla_specs,
              brms_specs = brms_specs
              ))
}


#' Helper function to predict fitted lme4/brms model
#' 
#' @description 
#' This function returns a df of cross-validated predictions
#' @details 
#' Expects numeric inputs
#' @param mod An [`g/lmerMod`] or [`brmsfit`] object.
#' @param pred_df_ls A list of 'test' data to make predictions for - each item in
#' list corresponds to an item in re_form.
#' @param re_form A [`list`] indicating the random effects variables to be
#' included in the cross-validation/insample predicitons.
#' @keywords model, evaluation, prediction
#' @examples 
#' \dontrun{
#'   preds <- predict_lme4_cv(mod_lme4_obj, pred_df_ls, list(NA, ~ (1|SS)))
#' } 
#' @noRd
predict_lme4_brms_cv <- function(mod, pred_df_ls, re_form){
  # Assertions
  assertthat::assert_that(
    inherits(mod, "lmerMod") || inherits(mod, "glmerMod") || inherits(mod, "brmsfit"),
    is.list(pred_df_ls),
    is.list(re_form),
    msg = "Input arguments are not of the correct type."
  )
  assertthat::assert_that(
    length(pred_df_ls) ==length(re_form),
    msg = "List of data for prediction has a different lenght to list of re_form options."
  )

  # predict
  for (i in seq_along(re_form)){
    # Check that each df has a column called "Random_effect_form"
    if (nrow(pred_df_ls[[i]]) > 0){
      assertthat::assert_that(
      "Random_effect_form" %in% colnames(pred_df_ls[[i]]),
      msg = "Random_effect_form column missing."
      )

      if (inherits(mod, "brmsfit")){
        tmp_pred <- stats::fitted(mod, pred_df_ls[[i]], 
                          scale = "response", # brms_specs$type, ??
                          re_form = re_form[[i]]) |>
            as.data.frame()
          pred_df_ls[[i]]$prediction <- tmp_pred$Estimate
          pred_df_ls[[i]]$prediction_q2.5 <- tmp_pred$Q2.5
          pred_df_ls[[i]]$prediction_q97.5 <- tmp_pred$Q97.5
      } else {
        # print(head(pred_df_ls[[i]]))
        # print(lme4::fixef((mod)))
        pred_df_ls[[i]]$prediction <- stats::predict(mod, pred_df_ls[[i]], 
                                              type = "response", # lme4_specs$type, ??
                                              re.form = re_form[[i]])
      }
    } 
  }

  out <- pred_df_ls |> dplyr::bind_rows()

  return(out)
}

