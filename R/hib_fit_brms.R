#' hib_fit_brms: Fit a brms style mixed-effects model
#' 
#' @description 
#' This function fits a brms model using the supplied database, formula and other
#' function arguments tuning the model fitting process.
#' 
#' @details 
#' This function takes the supplied database object (and various other arguments),
#' fits the specified model and returns a HibModelBRMS object.
#' @param obj A [`HibDatabase`] object or a list of them.
#' @param form A `formula` in the style of `brms` (as in `lme4`) - including fixed and random 
#' effects.
#' @param family A `family` function, e.g. gaussian, poisson - as allowed by brms
#' @param weight_specs A named `list` containing "weights" and "weight_type". \cr
#' - \code{weights}: a `numeric` vector containing the relative weights to
#' assign to each dataset/site, should sum to 1. \cr
#' - \code{weight_type}: a `character`, specifying either "dataset" or "site" - i.e. should
#' the weights argument assign relative weights to whole datasets, or for sites in
#' each dataset. \cr
#' 
#' Only applicable if a list of `HibDatabase`s are provided. \cr
#' The length of 'weights' should be equal to the number of `HibDatabase`s. \cr
#' Note: if using this argument, the supplied formula should have a response of 
#' the form `y | weights(data_weight) ~`.
#' The data_weight variable will be generated within the function call.
#' @param scale_specs A named `list` containing "covars_to_scale" and "scale_func". 
#' - `covars_to_scale` is a `character` vector of numeric covariates that should be 
#' scaled prior to modelling. 
#' - `scale_func` is a function that will be used to scale the specified 
#' covariates - if `NULL`, `base::scale()` will be used.
#' @param relevel_specs A `vector` where each item is: 
#' \code{"<covar_to_relevel_name>" = "<base_level>"}.
#' @param ... Any other named arguments that can be passed to `brms::brm()`, e.g., prior, data2
#' 
#' @return `HibModelBRMS`
#' 
#' @references 
#' * Paul-Christian Bürkner (2017). brms: An R Package for Bayesian Multilevel 
#'   Models Using Stan. Journal of Statistical Software, 80(1), 1-28. 
#'   doi:10.18637/jss.v080.i01
#' * Paul-Christian Bürkner (2018). Advanced Bayesian Multilevel Modeling with 
#'   the R Package brms. The R Journal, 10(1), 395-411. doi:10.32614/RJ-2018-017
#' * Paul-Christian Bürkner (2021). Bayesian Item Response Modeling in R with 
#'   brms and Stan. Journal of Statistical Software, 100(5), 1-54. 
#'   doi:10.18637/jss.v100.i05
#' 
#' @keywords model, fitting
#' 
#' @seealso
#' [HibModelBRMS]
#' \code{\link[brms]{brm}}
#' 
#' @examples 
#' \dontrun{
#'   mod <- hib_fit_brms(obj, y ~ x + ( 1 | a ), ...)
#' } 
#' @export
hib_fit_brms <- function(obj, 
                         form,
                         family = stats::gaussian,
                         weight_specs = list(weights = NULL,
                                             weight_type = NULL), # list(weights = c(a, b), weight_type = "dataset") # or "site"
                         scale_specs = list(covars_to_scale = NULL,
                                            scale_func = NULL),
                         relevel_specs = c(NULL),
                         ...
                         ){
  
  if(!isNamespaceLoaded("brms")) { attachNamespace("brms");requireNamespace('brms') }
  
  # Check inputs
  assertthat::assert_that(
    inherits(obj, "HibDatabase") || ( is.list(obj) & all(unlist(lapply(obj, function(x){inherits(x, "HibDatabase")}))) ) , # obj is a hib database
    inherits(form, "formula") || inherits(form, "brmsformula") || inherits(form, "mvbrmsformula"),
    is.function(family) || is.character(family), 
    is.list(weight_specs),
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
    # %%%
    mod_df_ls <- weight_model_data(mod_df_ls,
                                   "brms",
                                   weight_specs$weights,
                                   weight_specs$weight_type)
    # %%%
    
    mod_df <- dplyr::bind_rows(mod_df_ls)
    # ensure that mod_df has >1 row of data!
    assertthat::assert_that(
      nrow(mod_df) > 1
    )
    n_spp_taxa <- purrr::map2(.x = obj, 
                         .y = mod_df_ls, 
                         .f = function(.x, .y, resp = mod_resp){
                           calc_N_spp(.x, .y, resp)
                         }) #|> unlist()
    n_spp = lapply(n_spp_taxa, function(x){x$N_spp})
    n_taxa = lapply(n_spp_taxa, function(x){x$N_taxa})

    # N sites
    n_sites <- lapply(mod_df_ls, 
                      function(x){calc_N_sites(x)})
  }
  if (inherits(obj, "HibDatabase")){
    mod_df <- get_model_data(obj,
                             mod_resp,
                             mod_covars)
    # ensure that mod_df has >1 row of data!
    assertthat::assert_that(
      nrow(mod_df) > 1
    )
    n_spp_taxa <- calc_N_spp(obj, mod_df, mod_resp)
    n_spp = n_spp_taxa$N_spp
    n_taxa = n_spp_taxa$N_taxa
    n_sites <- calc_N_sites(mod_df)
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
  
  # backend = "cmdstanr"  - user can add
  # TODO:
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
  mod_obj <- HibModelBRMS$new(
    db_obj = obj, 
    mod_df = mod_df,
    brms_obj = m0,
    fixef_df = fixef_df,
    metr_ls = metr_ls,
    n_spp = n_spp,
    n_taxa = n_taxa,
    n_sites = n_sites,
    extra_args = list(...)
  )
  
  return(mod_obj)
}

