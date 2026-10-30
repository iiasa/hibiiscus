#' hib_fit_lme4: Fit an lme4 style mixed-effects model
#' 
#' @description 
#' This function fits an lme4 model using the supplied database, formula and other
#' function arguments tuning the model fitting process.
#' 
#' @details 
#' This function takes the supplied database object (and various other arguments),
#' fits the specified model and returns a HibModelLME4 object.
#' 
#' @param obj A [`HibDatabase`] object, or a list of databases.
#' @param form Model formula in the style of lme4 - including fixed and random 
#' effects.
#' @param family A `family` function, e.g. gaussian, poisson - as allowed by lme4.
#' Dafault is 'gaussian'.
#' @param REML A `logical` - if using a gaussian family, should the model be fitted 
#' using restricted maximum likelihood (REML, T) or maximum likelihood (F) - as in `lme4::lmer()`.
#' @param calc_CIs An indication of whether fixed effect confidence intervals (95%) 
#' should be calculated or not, and if so how - one of NULL, "Wald", "profile", "boot".
#' @param optimise_ranefs Not currently applicable - future dev
#' @param model_dredge Not currently applicable - future dev
#' @param model_average Not currently applicable - future dev
#' @param weight_specs A named `list` containing "weights" and "weight_type". \cr
#' - \code{weights}: a `numeric` vector containing the relative weights to
#' assign to each dataset/site, should sum to 1. \cr
#' - \code{weight_type}: a `character`, specifying either "dataset" or "site" - i.e. should
#' the weights argument assign relative weights to whole datasets, or for sites in
#' each dataset. \cr
#' 
#' Only applicable if a list of `HibDatabase`s are provided. \cr
#' The length of 'weights' should be equal to the number of `HibDatabase`s.
#' @param scale_specs A named `list` containing "covars_to_scale" and "scale_func". 
#' - `covars_to_scale` is a `character` vector of numeric covariates that should be 
#' scaled prior to modelling. 
#' - `scale_func` is a function that will be used to scale the specified 
#' covariates - if `NULL`, `base::scale()` will be used.
#' @param relevel_specs A `vector` where each item is: 
#' \code{"<covar_to_relevel_name>" = "<base_level>"}.
#' @param ... Any other named arguments that can be passed to `lme4::lmer()`/`lme4::glmer()`, 
#' e.g., na.action
#' 
#' @return `HibModelLME4`
#' 
#' @references 
#' * Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear 
#' Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48. 
#' doi:10.18637/jss.v067.i01.
#' 
#' @keywords model, fitting
#' 
#' @seealso
#' [HibModelLME4]
#' \code{\link[lme4]{lmer}}
#' \code{\link[lme4]{glmer}}
#' 
#' @examples 
#' \dontrun{
#'   mod <- hib_fit_lme4(obj, y ~ x + ( 1 | a ), ...)
#' } 
#' 
#' @export
hib_fit_lme4 <- function(obj, 
                         form,
                         family = stats::gaussian,
                         REML = TRUE,
                         calc_CIs = NULL, #,
                         optimise_ranefs = FALSE,
                         model_dredge = FALSE,
                         model_average = FALSE,
                         weight_specs = list(weights = NULL,
                                             weight_type = NULL), # list(weights = c(a, b), weight_type = "dataset") # or "site"
                         scale_specs = list(covars_to_scale = NULL,
                                            scale_func = NULL),
                         relevel_specs = c(NULL),
                         ...
                         ){
  
  # TODO:
  # Add functionality for optimal random effects structure - log-liklihood tests
  # Add functionality for dredging/testing a range of model params...?
  # Add functionality to return averaged model output
  
  rlang::check_installed(c("lme4", "MuMIn"),
                         reason = "to fit lme4 models and calculate model statistics.")
  
  # Add assertions
  assertthat::assert_that(
    inherits(obj, "HibDatabase") || ( is.list(obj) & all(unlist(lapply(obj, function(x){inherits(x, "HibDatabase")}))) ), # obj is a hib database
    # form can be a formula,
    is.function(family) || is.character(family), # family is character or function
    is.logical(REML), # REML is bool
    is.null(calc_CIs) || is.character(calc_CIs),   # calc_CIs is null or chr
    is.list(scale_specs),
    is.vector(relevel_specs) || is.null(relevel_specs),
    is.list(weight_specs)
  )


  # Check that calc_CIS is one of :
  if (!is.null(calc_CIs)){
    calc_CIs <- match.arg(calc_CIs, c("Wald", "profile", "boot"), several.ok = F)  
  }
  
  message(crayon::green("\nCollating data for model fitting...\n"))
  # Collate data from database - based on columns specified in the formula
  mod_covars <- all.vars(form[[3]])
  mod_resp <- all.vars(form[[2]])

  if (is.list(obj)){
    mod_df_ls <- purrr::map(.x = obj,
                            .f = function(.x, 
                                          resp = mod_resp,
                                          covars = mod_covars){
                              get_model_data(.x, resp, covars)
                            })
    
    # %%%
    # call weight calc function
    mod_df_ls <- weight_model_data(mod_df_ls,
                                   "lme4",
                                   weight_specs$weights,
                                   weight_specs$weight_type)
    
    mod_df <- dplyr::bind_rows(mod_df_ls,
                              .id = "dataset_modelling_id")
    # ensure that mod_df has >1 row of data!
    assertthat::assert_that(
      nrow(mod_df) > 1
    )
    n_spp_taxa <- purrr::map2(.x = obj, 
                         .y = mod_df_ls, 
                         .f = function(.x, .y, resp = mod_resp){
                           calc_N_spp(.x, .y, resp)
                         }) 
    n_spp = lapply(n_spp_taxa, function(x){x$N_spp})
    n_taxa = lapply(n_spp_taxa, function(x){x$N_taxa})
    # no need to unlist
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
    # N sites
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
  # get_all_vars(test_formula)
  # rlang::f_lhs(test_formula)
  # rlang::f_rhs(test_formula)
  
  message(crayon::green("Fitting model (this may take some time) ...\n"))
  
  # Specify lmer v glmer
  # first, check if family is a function or character...
  if (is.character(family)){
    if (family == "gaussian"){
      if (inherits(obj, "HibDatabase")){
        m0 <- lme4::lmer(form, 
                         data = mod_df,
                         REML = REML,
                         ...)
      }
      if (is.list(obj)){
        m0 <- lme4::lmer(form, 
                         data = mod_df,
                         REML = REML,
                         weights = data_weight,
                         ...)
      }
    }
    if (family != "gaussian"){
      if (inherits(obj, "HibDatabase")){
        m0 <- lme4::glmer(form, 
                          family = family,
                          data = mod_df,
                         ...)
      }
      
      if (is.list(obj)){
        m0 <- lme4::glmer(form, 
                         family = family,
                         data = mod_df,
                         weights = data_weight,
                         ...)
      }
    }
  }
  if (is.function(family)){
    if (family()$family == "gaussian"){
      if (inherits(obj, "HibDatabase")){
        m0 <- lme4::lmer(form, 
                         data = mod_df,
                         REML = REML,
                         ...)
      }
      if (is.list(obj)){
        m0 <- lme4::lmer(form, 
                         data = mod_df,
                         REML = REML,
                         weights = data_weight,
                         ...)
      }
    }
    if (family()$family != "gaussian"){
      if (inherits(obj, "HibDatabase")){
        m0 <- lme4::glmer(form, 
                          family = family,
                          data = mod_df,
                         ...)  
      }
      
      if (is.list(obj)){
        m0 <- lme4::glmer(form, 
                          family = family,
                          data = mod_df,
                         weights = data_weight,
                         ...)
      }
    }
  }
  

  # TODO: Do some checks - model convergence/model fitted
  

  # Collate model coefs../metrics - AICc, Rsq ...
  message(crayon::green("Collating model details ...\n"))
  
  fixef_df <- as.data.frame(lme4::fixef(m0))
  colnames(fixef_df) <- "Coefficient_value"
  fixef_df <- tibble::rownames_to_column(fixef_df, "Coefficient_name")
  
  if (!is.null(calc_CIs)){
    message(crayon::green("\tGetting confidence intervals ..."))
    confint_df <- lme4::confint.merMod(m0, parm = "beta_", method = calc_CIs)
    confint_df <- as.data.frame(confint_df)
    colnames(confint_df) <- c("loCI", "hiCI")
    confint_df$Coefficient_name <- rownames(confint_df)
    fixef_df <- dplyr::left_join(fixef_df,
                                 confint_df,
                                 by = "Coefficient_name")
  }
  
  metr_vec <- c(MuMIn::r.squaredGLMM(m0)[1,], stats::AIC(m0), stats::BIC(m0), MuMIn::AICc(m0))
  names(metr_vec) <- c("R2m", "R2c", "AIC", "BIC", "AICc")
  
  ci_type_ <- calc_CIs
  if (is.null(calc_CIs)){
    ci_type_ <- NA
  }

  # instantiate an HibModelLME4 object
  mod_obj <- HibModelLME4$new(
    db_obj = obj, 
    mod_df = mod_df,
    lme4_obj = m0,
    fixef_df = fixef_df,
    ci_type = ci_type_,
    metr_vec = metr_vec,
    n_spp = n_spp,
    n_taxa = n_taxa,
    n_sites = n_sites
  )
  
  return(mod_obj)
  
  # Future dev
  # Optional... extra wrappers...
  # If optimise_ranefs == F, fit m with ML
  # If optimise_ranefs == T, fit m with REML to first compare ranefs
  #   Then fit optimal ranef mod with ML
  
  # If dredge... dredge and model average
  
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal helper functions for data preparation
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' get_model_data: Collate model data from a HibDatabase
#' 
#' @description 
#' This function collects the data needed for the model fitting process, based on 
#' the supplied response and covariates
#' 
#' @param obj A `HibDatabase` object.
#' @param mod_resp A `character` indicating the response variable.
#' @param mod_covars A `character` indicating the model covariates.
#' 
#' @return data.frame
#' 
#' @keywords data processing
#' 
#' @examples 
#' \dontrun{
#'   mod_df <- get_model_data(obj, "Total_abundance", c("LUI"))
#' } 
#' @noRd
get_model_data <- function(obj,
                           mod_resp,
                           mod_covars){
  
  assertthat::assert_that(
    inherits(obj, "HibDatabase"),
    is.character(mod_resp),
    is.character(mod_covars)
    )

  within_site_resp <- purrr::map(obj$diversity_within$Data,
                                 colnames) |>
    unlist() |>
    unique() 
  between_site_resp <- purrr::map(obj$diversity_between$Data,
                                  colnames) |>
    unlist() |>
    unique() 
  
  # Get data for within-site/alpha diversity models.
  if (mod_resp %in% within_site_resp){
    
    # If not all covars are in obj$sites, as in a possible phylogenetically corr
    # brms model, subset to those covars available in sites and throw a warning
    all_covars <- unique(c(mod_covars, "SSBS"))
    missing_covars <- all_covars[!all_covars %in% colnames(obj$sites)]
    if (length(missing_covars) > 0){
      message(crayon::yellow(paste0("Warning: covariates specified in the model ", 
                                    "formula are missing; ",
                                    paste(missing_covars, colapse = ", "),
                                    "\n")))
    }
    all_covars_sub <- all_covars[all_covars %in% colnames(obj$sites)]
    mod_df <- obj$sites[,all_covars_sub]

    # find relevant obj$diversity$Data[[i]]   
    div_data_idx <- which(unlist(lapply(obj$diversity_within$Data,
                                        function(x, resp){
                                          tmp <- colnames(x)
                                          # any(grepl(resp, tmp))
                                          any(grepl(paste0("^", resp, "$"), tmp))
                                        },
                                        resp = mod_resp)))

    assertthat::assert_that(
      length(div_data_idx) == 1
    )
    
    mod_df <- mod_df |>
      dplyr::left_join(obj$diversity_within$Data[[div_data_idx]][,c(mod_resp, "SSBS")],
                      by = "SSBS")
    
    # Get data for between-site/beta diversity models.
  } else if (mod_resp %in% between_site_resp){

    all_covars <- unique(c(mod_covars, c("SSBS1", "SSBS2") ))
    missing_covars <- all_covars[!all_covars %in% colnames(obj$site_contrasts)]

    # If Dataset_name in missing covars, ignore it as a missing covar, 
    # add it to site_contrasts 
    incl_dat_nm <- FALSE
    if ("Dataset_name" %in% all_covars){
      missing_covars <- missing_covars[missing_covars != "Dataset_name"]
      incl_dat_nm <- TRUE
    }

    if (length(missing_covars) > 0){
      message(crayon::yellow(paste0("Warning: covariates specified in the model ", 
                                    "formula are missing; ",
                                    paste(missing_covars, colapse = ", "),
                                    "\n")))
    }
    all_covars_sub <- all_covars[all_covars %in% colnames(obj$site_contrasts)]
    mod_df <- obj$site_contrasts[,all_covars_sub]

    if (incl_dat_nm == TRUE){
      message(crayon::green("Fetching 'Dataset_name' information from HibDatabase\n"))
      mod_df <- mod_df |>
        dplyr::left_join(obj$sites[, c("SSBS", "Dataset_name")], 
                         by = c("SSBS1" = "SSBS"))
    }
    
    # find relevant obj$diversity$Data[[]] - which one contains the response...
    div_data_idx <- which(unlist(lapply(obj$diversity_between$Data,
                                        function(x, resp){
                                          tmp <- colnames(x)
                                          any(grepl(paste0("^", resp, "$"), tmp))
                                        },
                                        resp = mod_resp)))
    assertthat::assert_that(
      length(div_data_idx) == 1
    )
    
    mod_df <- mod_df |>
      dplyr::left_join(obj$diversity_between$Data[[div_data_idx]][,c(mod_resp, "SSBS1", "SSBS2")],
                      by = c("SSBS1", "SSBS2"))
    
  } else {
    message(crayon::magenta("Error: response variable specified is not available.\n"))
    stop()
  }
  
  # Drop rows with any NAs
  mod_df <- stats::na.omit(mod_df)
  
  return(mod_df)
  
}





#' prep_model_data: Prepare data for modelling
#' 
#' @description 
#' This function prepares the data for the model fitting based on user supplied
#' flags indicating the scaling of continuous variables and/or the re-leveling
#' of categorical covariates
#' 
#' @param df A `data.frame` to be used in model fitting.
#' @param mod_covars A `character` vector indicating the model covariates.
#' @param covars_to_scale A `character` vector indicating the response variables to scale.
#' @param scale_func A `function` to scale specified continuous variables.
#' @param covars_to_relevel A `character` vector containing categorical variables
#' to relevel.
#' @param base_levels A `character` vector containing baseline levels - in same
#' order as the variables in covars_to_relevel.
#' 
#' @return data.frame
#' 
#' @keywords data processing
#' 
#' @examples 
#' \dontrun{
#'   mod_df <- prep_model_data(mod_df, 
#'                             c(""),
#'                             scale,
#'                             c(""),
#'                             c("")
#' } 
#' @noRd
prep_model_data <- function(df, 
                            mod_covars,
                            covars_to_scale, scale_func, 
                            covars_to_relevel, base_levels){
  assertthat::assert_that(
    is.data.frame(df),
    is.character(mod_covars),
    is.character(covars_to_scale) || is.null(covars_to_scale),
    is.function(scale_func) || is.null(scale_func),
    is.character(covars_to_relevel) || is.null(covars_to_relevel),
    is.character(base_levels) || is.null(base_levels)
  )

  # FOR SAFETY (esp brms), convert any chr columns to fct (mod_covars)
  chr_cols_idx <- purrr::map(mod_covars, 
                             function(x, tmp = df){
                               is.character(tmp[[x]])}) |> 
    unlist()
  chr_cols <- mod_covars[chr_cols_idx]
  # print(chr_cols)
  for (c in chr_cols){
    df[[c]] <- as.factor(df[[c]])
  }
  
  # Check that covars are present
  if (!is.null(covars_to_scale)){
    assertthat::assert_that(
      # Covars exist
      all(covars_to_scale %in% colnames(df)),
      # check that scale_covars are numeric
      all(unlist(lapply(covars_to_scale, function(x, df_ = df){
        is.numeric(df[[x]])
      })))
    )

    # - check scale_func is not null
    if (is.null(scale_func)){
      message(crayon::yellow("Warning: scale_func not provided, using base R scale()."))
      scale_func <- scale
    }
    # scale,
    for (c in covars_to_scale){
      df[c] <- scale_func(df[[c]])
    }
  }
  if (!is.null(covars_to_relevel)){
    assertthat::assert_that(
      # Covars exist
      all(covars_to_relevel %in% colnames(df)),
      length(covars_to_relevel) == length(base_levels),
      # check that covars_to_relevel are chr/factor
      all(unlist(lapply(covars_to_relevel, function(x, df_ = df){
        ( is.character(df[[x]]) || is.factor(df[[x]]) )
      }))),
      # check factors exist for assoc covar
      all(unlist(purrr::map2(covars_to_relevel,
                             base_levels,
                             function(x, y, df_ = df){
                               y %in% df_[[x]]
                             })))
    )
    
    # TODO: could split the above assertions to give more meaningful error messages

    # drop unused factor levels
    df <- droplevels(df)

    # relevel
    for (i in 1:length(covars_to_relevel)){
      if (is.character(df[[covars_to_relevel[i]]])){
        df[[covars_to_relevel[i]]] <- factor(df[[covars_to_relevel[i]]])
      }
      df[[covars_to_relevel[i]]] <- stats::relevel(df[[covars_to_relevel[i]]],
                                             base_levels[i])
    }
  }
  
  return(df)
  
}



#' weight_model_data: Calculate weights for use in integrated modelling
#' 
#' @description 
#' This function calculates and assigns weights to rows in the modelling data
#' based on the provided arguments
#' 
#' @param df_ls A list of `data.frames`, each associated with a biodiversity dataset.
#' @param model_spec A `character` indicating if the model type is 'lme4', or 'brms'.
#' @param weights A `numeric` vector indicating the 'weights' the user wants to 
#' link to each dataset.
#' @param weight_type A `character` specifying the form of weighting. Either 
#' 'dataset' or 'site'. 
#' 
#' @return list
#' 
#' @keywords data processing
#' 
#' @examples 
#' \dontrun{
#'   mod_ls <- weight_model_data(df_ls, weights = c(0.5,0.5), weight_type = "dataset") 
#' } 
#' @noRd
weight_model_data <- function(df_ls,
                              model_spec,
                              weights,
                              weight_type){

  # TODO: need to add and check INLA options

  assertthat::assert_that(
    is.list(df_ls),
    is.character(model_spec),
    is.numeric(weights) || is.null(weights),
    is.character(weight_type) || is.null(weight_type)
  )
  
  assertthat::assert_that(
    all(unlist(lapply(df_ls, function(x){
        is.data.frame(x)
      })))
  )

  if (!is.null(weights)){
    assertthat::assert_that(
      length(df_ls) == length(weights)
    )
  }
  
  # TODO: assert that model_spec is within the potential model options
  # "lme4", "brms", ?"inla"?
  # check weight_type - NULL, dataset, site
  if (is.null(weight_type)){
    message(crayon::yellow("Warning: weight_type missing, setting to 'dataset'."))
    weight_type <- "dataset"
  }
  weight_type <- match.arg(weight_type, c("dataset", "site"), several.ok = F)  
  
  # check weights - NULL, numeric,.
  if (is.null(weights)){
    weights = rep(1, length(df_ls))
    message(crayon::yellow(paste0("Warning: weights not specified, using equal weight per ",
                                  weight_type, ".")))
  }
  
  
  # Need to know model type... - 
  # LME4 - weights sum to 1
  # BRMS - weights sum to N data points
  
  # Ensure weights sum to 1 (rescale)
  if (sum(weights) != 1){
    weights <- weights/sum(weights)  
  }
  
  
  if (weight_type == "dataset"){
    for (i in 1:length(df_ls)){
      df_ls[[i]]$data_weight <- weights[i] / nrow(df_ls[[i]])
    }
  }
  
  if (weight_type == "site"){
    # Total weight is equal to nrows of data - as above... but this is divided between
    # datasets/sites 
    # tot_dat_n <- sum(unlist(lapply(df_ls, nrow)))
    # weights <- weights / tot_dat_n
    for (i in 1:length(df_ls)){
      df_ls[[i]]$data_weight <- weights[i]
    }
    
    tot_w <- sum(unlist(lapply(df_ls, function(x){sum(x$data_weight)})))
    
    for (i in 1:length(df_ls)){
      df_ls[[i]]$data_weight <- df_ls[[i]]$data_weight / tot_w
    }
  }
  
  # If brms, just multiply above by tot_N
  if (model_spec == "brms"){
    tot_dat_n <- sum(unlist(lapply(df_ls, nrow)))
    for (i in 1:length(df_ls)){
      df_ls[[i]]$data_weight <- df_ls[[i]]$data_weight * tot_dat_n
    }
  }
  

  return(df_ls)
  
}


#' calc_N_sites: Calculate number of sites contributing to model data.frame
#' 
#' @description
#' This function calculates the total number of unique sites in a model data.frame
#' 
#' @param df A `data.frame`.
#' 
#' @return numeric
#' 
#' @keywords data processing
#' 
#' @examples
#' \dontrun{
#'   N_site <- calc_N_sites(df)
#' }
#' @noRd
calc_N_sites <- function(df){
  if ("SSBS" %in% colnames(df)){
    n_sites <- length(unique(df$SSBS))
  } else if ("SSBS1" %in% colnames(df)){
    n_sites <- length(unique(c(df$SSBS1, df$SSBS2)))
  }
  return(n_sites)                
}

