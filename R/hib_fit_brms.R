#' Fit a brms style mixed-effects model
#' 
#' @description 
#' This function fits a brms model using the supplied database, formula and other
#' function arguments tuning the model fitting process
#' @details 
#' Expects a database object...
#' @param obj A [`Database`] object or a list of them.
#' @param form A [`formula`] in the style of brms (as in lme4) - including fixed and random 
#' effects.
#' @param family A [`family`] function, e.g. gaussian, poisson - as allowed by brms
#' @param weight_specs A [`list`] containing "weights" and "weight_type".
#' weights is [`numeric`] vector containing the relative weights to
#' assign to each dataset/site, should sum to 1. 
#' weight_type is a character, specifying either "dataset" or "site" - i.e. should
#' the weights argument assign relative weights to whole datasets, or for sites in
#' each dataset.
#' Note: if using this argument, the supplied formula should have a response of 
#' the form [`y | weights(data_weight) ~`].
#' The data_weight variable will be generated within the function call.
#' @param scale_specs A [`list`] containing "covars_to_scale" and "scale_func". 
#' covars_to_scale is a [`character`] vector of numeric covariates that should be 
#' scaled prior to modelling. scale_func is a function that will be used to scale 
#' the specified covariates - if NULL, base::scale() will be used.
#' @param relevel_specs A named [`vector`] where each item is:
#' "<covar_to_relevel_name>" = "<base_level>".
#' @param ... Any other named arguments that can be passed to brms::brm(), e.g. prior, data2
#' @keywords model, fitting
#' @examples 
#' \dontrun{
#'   mod <- hib_fit_brms(obj, y ~ x + ( 1 | a ), ...)
#' } 
#' @export
hib_fit_brms <- function(obj, 
                         form,
                         family = gaussian,
                         # dataset_weights = NULL,
                         weight_specs = list(weights = NULL,
                                             weight_type = NULL), # list(weights = c(a, b), weight_type = "dataset") # or "site"
                         scale_specs = list(covars_to_scale = NULL,
                                            scale_func = NULL),
                         relevel_specs = c(NULL),
                         ...
                         ){
  
  if(!isNamespaceLoaded("brms")) { attachNamespace("brms");requireNamespace('brms') }
  
  # Add assertions
  assertthat::assert_that(
    inherits(obj, "Database") || ( is.list(obj) & all(unlist(lapply(obj, function(x){inherits(x, "Database")}))) ) , # obj is a hib database
    # form can be a formula,
    inherits(form, "formula") || inherits(form, "brmsformula") || inherits(form, "mvbrmsformula"),
    is.function(family) || is.character(family), # family is character or function
    # is.numeric(dataset_weights) || is.null(dataset_weights),
    is.list(weight_specs), #|| is.null(weight_specs)
    is.list(scale_specs),
    is.vector(relevel_specs) || is.null(relevel_specs)
  )
  
  message(crayon::green("\nCollating data for model fitting...\n"))
  
  
  # Collate data from database - based on columns specified in the formula
  if (inherits(form, "formula")){
    mod_covars <- all.vars(form[[3]])
    mod_resp <- all.vars(form[[2]])
  }
  
  
  if ((inherits(form, "brmsformula") | inherits(form, "mvbrmsformula"))){
    # stop("use of mv/brmsformulas is in progress!!")
    # First, try to use brmsterms and allterms output
    tmp <- try(brms::brmsterms(form, resp_rhs_all = F))
    # If that breaks, just use the first formula, and get terms as above...
    if (inherits(tmp, "try-error")){
      tmp_form <- form$formula
      message(crayon::yellow("Warning: terms from brmsformula were not peoperly parsed.. \n"))
    } else {
      tmp_form <- tmp$allvars
    }
    
    mod_covars <- all.vars(tmp_form[[3]])
    mod_resp <- all.vars(tmp_form[[2]])
  }
  
  # Drop data_weights from response
  mod_resp <- mod_resp[mod_resp != "data_weight"]
  
  # %%%%%
  # if weights is not null, ensure sum to 1 (just rescale)
  # if (!is.null(dataset_weights)){
  #   assertthat::assert_that(
  #     length(dataset_weights) == length(obj)
  #   )
  #   if (sum(dataset_weights) != 1){
  #     dataset_weights <- dataset_weights/sum(dataset_weights)  
  #   }
  # }
  # # if obj is a list, and weights is null, set to equal weight per dataset
  # if (is.list(obj) & is.null(dataset_weights)){
  #   dataset_weights <- rep(1/length(obj), length(obj)) 
  # }
  # 
  # %%%%%
  # Check weights - could move to weights func
  # if (is.list(obj) & !is.null(weight_specs$weights)){
  #   # Check weight_type is not null
  #   if (is.null(weight_specs$weight_type)){
  #     message(crayon::yellow("Warning: weight_type missing, setting to 'dataset'."))
  #     weight_specs$weight_type <- "dataset"
  #   }
  # }
  
  # if (is.list(obj) & is.null(weight_specs$weights)){
  #   # Default option of equal weight - if weights is null
  #   weight_specs <- list(weights = rep(1, length(obj)), weight_type = weight_specs$weight_type)
    
  #   # Check weight_type 
  #   if (!is.null(weight_specs$weight_type)){
  #     message(crayon::yellow(paste0("Warning: weights not specified, using equal weight per ",
  #                                   weight_specs$weight_type, ".")))
  #   }
  #   if (is.null(weight_specs$weight_type)){
  #     message(crayon::yellow("Warning: weights not specified and weight_type missing, using equal weight per dataset."))
  #     weight_specs$weight_type <- "dataset"
  #   }
  # }

  # %%%%%
  
  # %%%%%
  # TODO:
  # Check for covariance structures in model formula (mod_covars) and data2
  # if data2 exists, check the names
  # if names of data2 match something in mod_covars,  drop that value from mod_covars
  if ("data2" %in% names(list(...))){
    tmp_dots_ <- list(...)
    dat2_nms <- names(dots_[["data2"]])
    mod_covars <- mod_covars[!mod_covars %in% dat2_nms]
    rm(tmp_dots_)
  }
  # %%%%%
  
  
  # TODO:
  # What about a formula of the type 
  # y ~ 0 + Intercept + b
  # Intercept is an internally generated brms thing...
  if ("Intercept" %in% mod_covars){
    mod_covars <- mod_covars[!mod_covars == "Intercept"]
  }
  
  # %%%%%
  if (is.list(obj)){
    mod_df_ls <- purrr::map(.x = obj,
                            .f = function(.x, 
                                          resp = mod_resp,
                                          covars = mod_covars){
                              get_model_data(.x, resp, covars)
                            })
    
    # add weights
    # for (i in 1:length(obj)){
    #   mod_df_ls[[i]]$data_weight <- dataset_weights[i] / nrow(mod_df_ls[[i]])
    #   # This seems wrong... - 
    #   # For brms, weight is taken literally, minimum should be 1... maybe
    #   # 3 (initial) possible approaches to deriving data-point weights
    #   # - as above, sum to 1 over all datapoints
    #   # - sum to total number of datapoints, assign p1 * N to data1 and p2 * N to data2
    #   #   then divide these amounts by N_data_i - THIS SEEMS MOST APPROPRIATE
    #   # - set baseline weight to 1 for lower weighted dataset..(or higher weighted dataset)
    #   #   then adjust weights of other dataset/s accordingly...
    # }
    
    # %%%
    # n_df_vec <- lapply(mod_df_ls, nrow) |> unlist()
    # 
    # prop_weight_vec <- dataset_weights * sum(n_df_vec)
    # df_weight_vec <- prop_weight_vec / n_df_vec
    # 
    # for (i in 1:length(obj)){
    #   mod_df_ls[[i]]$data_weight <- df_weight_vec[i]
    # }
    
    # %%%
    mod_df_ls <- weight_model_data(mod_df_ls,
                                   "brms",
                                   weight_specs$weights,
                                   weight_specs$weight_type)
    # %%%
    
    mod_df <- dplyr::bind_rows(mod_df_ls)
    
    n_spp <- purrr::map2(.x = obj, 
                         .y = mod_df_ls, 
                         .f = function(.x, .y, resp = mod_resp){
                           calc_N_spp(.x, .y, resp)
                         }) |> unlist()
  }
  if (inherits(obj, "Database")){
    mod_df <- get_model_data(obj,
                             mod_resp,
                             mod_covars)
    n_spp <- calc_N_spp(obj, mod_df, mod_resp)
  }
  
  # %%%
  mod_df <- prep_model_data(mod_df,
                            mod_covars,
                            scale_specs$covars_to_scale,
                            scale_specs$scale_func,
                            names(relevel_specs),
                            unname(relevel_specs))
  # %%%

  message(crayon::green("Fitting model (this may take some time) ..."))
  
  m0 <- brms::brm(form, 
                  data = mod_df,
                  family = family,
                  ...)
  
  # Do some checks - model convergence/model fitted
  # Rhat...
  

  # Collate model coefs../metrics - AICc, Rsq ...
  message(crayon::green("Collating model details ..."))
  
  fixef_df <- as.data.frame(brms::fixef(m0))
  colnames(fixef_df) <- c("Coefficient_value", "Coefficient_error", "loCI", "hiCI")
  fixef_df <- tibble::rownames_to_column(fixef_df, "Coefficient_name")
  
  # waic, loo# Get loo and waic
  loo_metr <- brms::loo(m0)
  waic_metr <- brms::waic(m0)
  
  loo_metr$estimates
  waic_metr$estimates
  
  metr_ls <- list("loo" = loo_metr,
                  "waic" = waic_metr)

  # instantiate a brms_mod class
  mod_obj <- hib_brms_model$new(
    db_obj = obj, 
    mod_df = mod_df,
    brms_obj = m0,
    fixef_df = fixef_df,
    metr_ls = metr_ls,
    n_spp = n_spp,
    extra_args = list(...)
  )
  
  return(mod_obj)
}