#' Fit an lme4 style mixed-effects model
#' 
#' @description 
#' This function fits an lme4 model using the supplied database, formula and other
#' function arguments tuning the model fitting process
#' @details 
#' Expects a database object...
#' @param obj A [`Database`] object, or a list of databases.
#' @param form A [`formula`] in the style of lme4 - including fixed and random 
#' effects.
#' @param family A [`family`] function, e.g. gaussian, poisson - as allowed by lme4
#' @param REML A [`boolean`] - if using a gaussian family, should model be fitted 
#' useing REML or ML - as in lme4
#' @param calc_CIs An indication of whether fixed effect confidence intervals (95%) 
#' should be calculated or not, and if so how - one of NULL, "Wald", "profile", "boot"
#' @param optimise_ranefs Not currently applicable - future dev
#' @param model_dredge Not currently applicable - future dev
#' @param model_average Not currently applicable - future dev
#' @param dataset_weights A [`numeric`] vector containing the relative weights to
#' assign to each dataset, should sum to 1. Note: if using this argument, the 
#' supplied formula should have a response of the form [`y | weights(data_weight) ~`]
#' The data_weight variable will be generated within the function call so as to 
#' assign each observation per dataset an equal fraction of that dataset's weight.
#' @keywords model, fitting
#' @examples 
#' \dontrun{
#'   mod <- hib_fit_lme4(obj, y ~ x + ( 1 | a ), ...)
#' } 
#' @export
hib_fit_lme4 <- function(obj, 
                         form,
                         family = gaussian,
                         REML = T,
                         calc_CIs = NULL, #,
                         optimise_ranefs = F,
                         model_dredge = F,
                         model_average = F,
                         dataset_weights = NULL
                         ){
  
  # TODO:
  # Add functionality for optimal random effects structure - log-liklihood tests
  # Add functionality for dredging/testing a range of model params...?
  # Add functionality to return averaged model output
  
  # Check whether lme4 package is available
  # TODO: figure out check_package...
  # check_package('lme4')
  if(!isNamespaceLoaded("lme4")) { attachNamespace("lme4");requireNamespace('lme4') }
  
  # file.path("lme4")
  # rcmdcheck::rcmdcheck(path = "lme4",
  #                      libpath = .libPaths())
  
  # Add assertions
  assertthat::assert_that(
    inherits(obj, "Database") || is.list(obj), # obj is a hib database
    # form can be a formula,
    is.function(family) || is.character(family), # family is character or function
    is.logical(REML), # REML is bool
    is.null(calc_CIs) || is.character(calc_CIs),   # calc_CIs is null or chr
    is.null(dataset_weights) || is.numeric(dataset_weights)
  )


  # Check that calc_CIS is one of :
  if (!is.null(calc_CIs)){
    calc_CIs <- match.arg(calc_CIs, c("Wald", "profile", "boot"), several.ok = F)  
  }
  
  message(crayon::green("Collating data for model fitting...\n"))
  # Collate data from database - based on columns specified in the formula
  
  # %%%%%
  # if weights is not null, ensure sum to 1 (just rescale)
  if (!is.null(dataset_weights)){
    if (sum(dataset_weights) != 1){
      dataset_weights <- dataset_weights/sum(dataset_weights)  
    }
  }
  # if obj is a list, and weights is null, set to equal weight per dataset
  if (is.list(obj) & is.null(dataset_weights)){
    dataset_weights <- rep(1/length(obj), length(obj)) 
  }
  # %%%%%
  
  mod_covars <- all.vars(form[[3]])
  mod_resp <- all.vars(form[[2]])
  
  
  if (is.list(obj)){
    mod_df_ls <- purrr::map(.x = obj,
                            .f = function(.x, 
                                          resp = mod_resp,
                                          covars = mod_covars){
                              get_model_data(.x, resp, covars)
                            })
    
    # add weights
    for (i in 1:length(obj)){
      mod_df_ls[[i]]$data_weight <- dataset_weights[i] / nrow(mod_df_ls[[i]])
    }
    mod_df <- dplyr::bind_rows(mod_df_ls,
                              .id = "dataset_modelling_id")
    
    n_spp <- purrr::map2(.x = obj, 
                         .y = mod_df_ls, 
                         .f = function(.x, .y, resp = mod_resp){
                           calc_N_spp(.x, .y, resp)
                         }) %>% unlist()
  }
  if (inherits(obj, "Database")){
    mod_df <- get_model_data(obj,
                             mod_resp,
                             mod_covars)
    n_spp <- calc_N_spp(obj, mod_df, mod_resp)
  }

  
  # TODO:
  # Re-level LU/UI/LUI columns - optional... what if they include cannot decide... 
  # - or have other, groupings - primary, secondary...
  # lu_levels_vec <- c(
  #   "Primary vegetation",
  #   "Primary forest",
  #   "Primary non-forest",
  #   "Young secondary vegetation",
  #   "Intermediate secondary vegetation",
  #   "Mature secondary vegetation",
  #   "Secondary vegetation (indeterminate age)", 
  #   "Plantation forest",
  #   "Pasture",
  #   "Cropland",
  #   "Urban" ,
  #   "Cannot decide"
  # )
  # ui_levels_vec <- c(
  #   "Minimal use",
  #   "Light use",
  #   "Intense use",
  #   "Cannot decide"
  # )
  # lui_levels_vec <- c(
  #   "Primary vegetation_Minimal use",
  #   "Primary vegetation_Light use",
  #   "Primary vegetation_Intense use",
  #   
  #   "Primary forest_Minimal use",
  #   "Primary forest_Light use",
  #   "Primary forest_Intense use",
  #   
  #   "Primary non-forest_Minimal use",
  #   "Primary non-forest_Light use",
  #   "Primary non-forest_Intense use",
  #   
  #   "Mature secondary vegetation_Minimal use",
  #   "Mature secondary vegetation_Light use",
  #   "Mature secondary vegetation_Intense use",  
  #   
  #   "Intermediate secondary vegetation_Minimal use",
  #   "Intermediate secondary vegetation_Light use",
  #   "Intermediate secondary vegetation_Intense use",
  #   
  #   "Young secondary vegetation_Minimal use",
  #   "Young secondary vegetation_Light use",
  #   "Young secondary vegetation_Intense use",
  #   
  #   "Secondary vegetation (indeterminate age)_Minimal use",
  #   "Secondary vegetation (indeterminate age)_Light use",
  #   "Secondary vegetation (indeterminate age)_Intense use",
  #   
  #   "Plantation forest_Minimal use",
  #   "Plantation forest_Light use",
  #   "Plantation forest_Intense use",
  #   
  #   "Pasture_Minimal use",
  #   "Pasture_Light use",
  #   "Pasture_Intense use",
  #   
  #   "Cropland_Minimal use", 
  #   "Cropland_Light use",
  #   "Cropland_Intense use",
  #   
  #   "Urban_Minimal use",
  #   "Urban_Light use",                                                     
  #   "Urban_Intense use" 
  # )
  
  
  
  # get_all_vars(test_formula)
  # rlang::f_lhs(test_formula)
  # rlang::f_rhs(test_formula)
  
  # Specify lmer v glmer
  message(crayon::green("Fitting model (this may take some time) ..."))
  
  
  # first, check if family is a function or character...
  if (is.character(family)){
    if (family == "gaussian"){
      if (inherits(obj, "Database")){
        m0 <- lme4::lmer(form, 
                         data = mod_df,
                         REML = REML)
      }
      if (is.list(obj)){
        m0 <- lme4::lmer(form, 
                         data = mod_df,
                         REML = REML,
                         weights = data_weight)
      }
    }
    if (family != "gaussian"){
      if (inherits(obj, "Database")){
        m0 <- lme4::glmer(form, 
                          family = family,
                          data = mod_df)
      }
      
      if (is.list(obj)){
        m0 <- lme4::glmer(form, 
                         family = family,
                         data = mod_df,
                         weights = data_weight)
      }
    }
  }
  if (is.function(family)){
    if (family()$family == "gaussian"){
      if (inherits(obj, "Database")){
        m0 <- lme4::lmer(form, 
                         data = mod_df,
                         REML = REML)
      }
      if (is.list(obj)){
        m0 <- lme4::lmer(form, 
                         data = mod_df,
                         REML = REML,
                         weights = data_weight)
      }
    }
    if (family()$family != "gaussian"){
      if (inherits(obj, "Database")){
        m0 <- lme4::glmer(form, 
                          family = family,
                          data = mod_df)  
      }
      
      if (is.list(obj)){
        m0 <- lme4::glmer(form, 
                          family = family,
                          data = mod_df,
                         weights = data_weight)
      }
    }
  }
  

  # Do some checks - model convergence/model fitted
  
  # class(m0)
  # Collate model coefs../metrics - AICc, Rsq ...
  message(crayon::green("Collating model details ..."))
  
  fixef_df <- as.data.frame(lme4::fixef(m0))
  colnames(fixef_df) <- "Coefficient_value"
  fixef_df <- tibble::rownames_to_column(fixef_df, "Coefficient_name")
  # print(calc_CIs)
  if (!is.null(calc_CIs)){
    message(crayon::green("Getting confidence intervals ..."))
    confint_df <- lme4::confint.merMod(m0, parm = "beta_", method = calc_CIs)
    confint_df <- as.data.frame(confint_df)
    colnames(confint_df) <- c("loCI", "hiCI")
    confint_df$Coefficient_name <- rownames(confint_df)
    fixef_df <- dplyr::left_join(fixef_df,
                                 confint_df)
  }
  
  metr_vec <- c(MuMIn::r.squaredGLMM(m0)[1,], AIC(m0), BIC(m0), MuMIn::AICc(m0))
  names(metr_vec) <- c("R2m", "R2c", "AIC", "BIC", "AICc")
  
  ci_type_ <- calc_CIs
  if (is.null(calc_CIs)){
    ci_type_ <- NA
  }
  # instantiate an lme4_mod class
  mod_obj <- hib_lme4_model$new(
    db_obj = obj, 
    mod_df = mod_df,
    lme4_obj = m0,
    fixef_df = fixef_df,
    ci_type = ci_type_,
    metr_vec = metr_vec,
    n_spp = n_spp
  )
  
  return(mod_obj)
  # 
  
  # Optional... extra wrappers...
  # If optimise_ranefs == F, fit m with ML
  # If optimise_ranefs == T, fit m with REML to first compare ranefs
  #   Then fit optimal ranef mod with ML
  
  # If dredge... dredge and model average
  
}


#' Collate model data for mod_df
#' 
#' @description 
#' This function collest the data needed for the model fitting process, based on 
#' the supplied response and covariables
#' @details 
#' Expects a database object...
#' @param obj A [`Database`] object.
#' @param mod_resp A [`character`] indicating the response variable.
#' @param mod_covars A [`character`] indicating the model covariates.
#' @keywords data processing
#' @examples 
#' \dontrun{
#'   mod_df <- get_model_data(obj, "Total_abundance", c("LUI"))
#' } 
#' @noRd
get_model_data <- function(obj,
                           mod_resp,
                           mod_covars){
  
  within_site_resp <- purrr::map(obj$diversity_within$Data,
                                 colnames) %>%
    unlist() %>%
    unique() 
  between_site_resp <- purrr::map(obj$diversity_between$Data,
                                  colnames) %>%
    unlist() %>%
    unique() 
  
  if (mod_resp %in% within_site_resp){
    # If not all covars are in obj$sites, as in a possible phylogenetically corr
    # brms model, this will fail!!
    # Could subset to those covars available in sites - as in INLA
    # And throw a warning?!
    
    mod_df <- obj$sites[,unique(c(mod_covars, c("SS", "SSB", "SSBS")))]
    # find relevant obj$diversity$Data[[i]] 
    div_data_idx <- which(unlist(lapply(obj$diversity_within$Data,
                                        function(x, resp){
                                          tmp <- colnames(x)
                                          any(grepl(resp, tmp))
                                        },
                                        resp = mod_resp)))
    assertthat::assert_that(
      length(div_data_idx) == 1
    )
    
    mod_df <- mod_df %>%
      dplyr::left_join(obj$diversity_within$Data[[div_data_idx]][,c(mod_resp, "SSBS")])
    
  } else if (mod_resp %in% between_site_resp){
    mod_df <- obj$site_contrasts[,unique(c(mod_covars, c("SSBS1", "SSBS2", "SS", "SSB")))]
    
    # find relevant obj$diversity$Data[[]] - which one contains the response...
    div_data_idx <- which(unlist(lapply(obj$diversity_between$Data,
                                        function(x, resp){
                                          tmp <- colnames(x)
                                          any(grepl(resp, tmp))
                                        },
                                        resp = mod_resp)))
    assertthat::assert_that(
      length(div_data_idx) == 1
    )
    
    mod_df <- mod_df %>%
      dplyr::left_join(obj$diversity_between$Data[[div_data_idx]][,c(mod_resp, "SSBS1", "SSBS2")])
    
  } else {
    message(crayon::magenta("Error: response variable specified is not available.\n"))
    stop()
  }
  
  mod_df <- na.omit(mod_df)
  
  return(mod_df)
  
}