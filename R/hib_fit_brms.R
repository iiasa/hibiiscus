#' Fit a brms style mixed-effects model
#' 
#' @description 
#' This function fits a brms model using the supplied database, formula and other
#' function arguments tuning the model fitting process
#' @details 
#' Expects a database object...
#' @param obj A [`Database`] object.
#' @param form A [`formula`] in the style of brms (as in lme4) - including fixed and random 
#' effects.
#' @param family A [`family`] function, e.g. gaussian, poisson - as allowed by brms
#' @keywords model, fitting
#' @examples 
#' \dontrun{
#'   mod <- hib_fit_brms(obj, y ~ x + ( 1 | a ), ...)
#' } 
#' @export
hib_fit_brms <- function(obj, 
                         form,
                         family = gaussian
                         # TODO: 
                         # - add a list of model fit tuning params, named as in brms
                         # - data2 - for autocorrelation... - 
                         # PASS ALL EXTRA ARGS AS "..."?!
                         # GATHER IN FUNCTION
                         ){
  
  if(!isNamespaceLoaded("brms")) { attachNamespace("brms");requireNamespace('brms') }
  
  # Add assertions
  assertthat::assert_that(
    inherits(obj, "Database"), # obj is a hib database
    # form can be a formula,
    inherits(form, "formula"),
    is.function(family) || is.character(family) # family is character or function
  )
  
  message(crayon::green("Collating data for model fitting...\n"))
  
  
  # Collate data from database - based on columns specified in the formula
  
  
  mod_covars <- all.vars(form[[3]])
  mod_resp <- all.vars(form[[2]])
  

  # consider diversity_within and diversity_between
  # check where response variable is - diversity_within v diversity_between
  # if response is in diversity_within, need to use sites df for covars
  # if response is in diversity_between, need to use site_contrasts df for covars
  
  within_site_resp <- purrr::map(obj$diversity_within$Data,
                                 colnames) %>%
    unlist() %>%
    unique() 
  between_site_resp <- purrr::map(obj$diversity_between$Data,
                                  colnames) %>%
    unlist() %>%
    unique() 
  
  if (mod_resp %in% within_site_resp){
    # print("getting model data")
    mod_df <- obj$sites[,unique(c(mod_covars, c("SS", "SSB", "SSBS")))]
    # find relevant obj$diversity$Data[[]] - which one contains the response...
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
    
    mod_df <- na.omit(mod_df)
    
    n_spp <- obj$diversity_within$N_species_total[[div_data_idx]]
    
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
    
    mod_df <- na.omit(mod_df)
    
    n_spp <- obj$diversity_between$N_species_total[[div_data_idx]]
    
  } else {
    message(crayon::magenta("Error: response variable specified is not available.\n"))
    stop()
  }
  
  
  
  message(crayon::green("Fitting model (this may take some time) ..."))
  
  m0 <- brms::brm(form, 
                  data = mod_df,
                  family = family)
  
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