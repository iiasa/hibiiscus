#' Fit an lme4 stype mixed-effects model
#' 
#' @description 
#' This function fits an lme4 model using the supplied database, formula and other
#' function arguments tuning the model fitting process
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @param form A [`formula`] in the style of lme4 - including fixed and random effects.
#' 
#' @keywords model, fitting
#' @examples 
#' \dontrun{
#'   mod <- hib_fit_lme4(obj, ...)
#' } 
#' @export
hib_fit_lme4 <- function(obj, 
                         form,
                         family = gaussian,
                         REML = T,
                         # link = "identity",
                         optimise_ranefs = F,
                         model_dredge = F,
                         model_average = F){
  
  # Check whether lme4 package is available
  # TODO: figure out check_package...
  # check_package('lme4')
  if(!isNamespaceLoaded("lme4")) { attachNamespace("lme4");requireNamespace('lme4') }
  
  # file.path("lme4")
  # rcmdcheck::rcmdcheck(path = "lme4",
  #                      libpath = .libPaths())
  
  # TODO: add progress messages
  message(crayon::green("Collating data for model fitting...\n"))
  # Collate data from database - based on columns specified in the formula
  mod_covars <- all.vars(form[[3]])
  mod_resp <- all.vars(form[[2]])
  
  mod_df <- obj$sites[,unique(c(mod_covars, c("SS", "SSB", "SSBS")))]
  
  # find relevant obj$diversity$Data[[]] - which one contains the response...
  div_data_idx <- which(unlist(lapply(obj$diversity$Data,
                                      function(x, resp){
                                        tmp <- colnames(x)
                                        any(grepl(resp, tmp))
                                      },
                                      resp = mod_resp)))
  assertthat::assert_that(
    length(div_data_idx) == 1
  )
  
  mod_df <- mod_df %>%
    dplyr::left_join(obj$diversity$Data[[div_data_idx]][,c(mod_resp, "SSBS")])
  
  mod_df <- na.omit(mod_df)
  
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
      m0 <- lme4::lmer(form, 
                       data = mod_df,
                       REML = REML)
    }
    if (family != "gaussian"){
      m0 <- lme4::glmer(form, 
                        family = family,
                        data = mod_df)
    }
  }
  if (is.function(family)){
    if (family()$family == "gaussian"){
      m0 <- lme4::lmer(form, 
                       data = mod_df,
                       REML = REML)
    }
    if (family()$family != "gaussian"){
      m0 <- lme4::glmer(form, 
                        family = family,
                        data = mod_df)
    }
  }
  

  # Do some checks - model convergence/model fitted
  
  # class(m0)
  # Collate model coefs../metrics - AICc, Rsq ...
  message(crayon::green("Collating model details ..."))
  
  fixef_df <- as.data.frame(lme4::fixef(m0))
  colnames(fixef_df) <- "Coefficient_value"
  fixef_df <- tibble::rownames_to_column(fixef_df, "Coefficient_name")
  
  # confint(m0, method = "Wald")
  # confint(m0, method="boot", nsim=1000, parm="beta_")
  # lme4::confint.merMod(m0)
  
  metr_vec <- c(MuMIn::r.squaredGLMM(m0)[1,], summary(m0)$AICtab[c("AIC", "BIC")], MuMIn::AICc(m0))
  names(metr_vec) <- c("R2m", "R2c", "AIC", "BIC", "AICc")
  
  n_spp <- obj$diversity$N_species_total[[div_data_idx]]
  
  # instantiate a lme4_mod class
  mod_obj <- hib_lme4_model$new(
    db_obj = obj, 
    mod_df = mod_df,
    lme4_obj = m0,
    fixef_df = fixef_df,
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
  
  # 
  
  # %%%%%
  # testing
  # form <- sqrt(Rescaled_total_abundance) ~ LUI + (1 | SS) + (1 | SSB)
  # m0 <- lme4::lmer(form, 
  #                  data = mod_df,
  #                  REML = F)
  # %%%%%

  
}