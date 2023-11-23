#' Fit an lme4 style mixed-effects model
#' 
#' @description 
#' This function fits an lme4 model using the supplied database, formula and other
#' function arguments tuning the model fitting process
#' @details 
#' Expects a database object and various other arguments.
#' @param obj A [`Database`] object, or a list of databases.
#' @param form A [`formula`] in the style of lme4 - including fixed and random 
#' effects.
#' @param family A [`family`] function, e.g. gaussian, poisson - as allowed by lme4
#' @param REML A [`logical`] - if using a gaussian family, should model be fitted 
#' using REML (T) or ML (F) - as in lme4
#' @param calc_CIs An indication of whether fixed effect confidence intervals (95%) 
#' should be calculated or not, and if so how - one of NULL, "Wald", "profile", "boot"
#' @param optimise_ranefs Not currently applicable - future dev
#' @param model_dredge Not currently applicable - future dev
#' @param model_average Not currently applicable - future dev
#' @param weight_specs A [`list`] containing "weights" and "weight_type".
#' weights is a [`numeric`] vector containing the relative weights to
#' assign to each dataset/site, should sum to 1. 
#' weight_type is a [`character`], specifying either "dataset" or "site" - i.e. should
#' the weights argument assign relative weights to whole datasets, or for sites in
#' each dataset.
#' @param scale_specs A [`list`] containing "covars_to_scale" and "scale_func". 
#' covars_to_scale is a [`character`] vector of numeric covariates that should be 
#' scaled prior to modelling. 
#' scale_func is a function that will be used to scale 
#' the specified covariates - if NULL, base::scale() will be used.
#' @param relevel_specs A named [`vector`] where each item is:
#' "<covar_to_relevel_name>" = "<base_level>".
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
                         # dataset_weights = NULL,
                         weight_specs = list(weights = NULL,
                                             weight_type = NULL), # list(weights = c(a, b), weight_type = "dataset") # or "site"
                         scale_specs = list(covars_to_scale = NULL,
                                            scale_func = NULL),
                         relevel_specs = c(NULL)
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
    inherits(obj, "Database") || ( is.list(obj) & all(unlist(lapply(obj, function(x){inherits(x, "Database")}))) ), # obj is a hib database
    # form can be a formula,
    is.function(family) || is.character(family), # family is character or function
    is.logical(REML), # REML is bool
    is.null(calc_CIs) || is.character(calc_CIs),   # calc_CIs is null or chr
    # is.null(dataset_weights) || is.numeric(dataset_weights),
    is.list(scale_specs),
    is.vector(relevel_specs) || is.null(relevel_specs),
    is.list(weight_specs) #|| is.null(weight_specs)
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
    # TODO : call weight calc function
    mod_df_ls <- weight_model_data(mod_df_ls,
                                   "lme4",
                                   weight_specs$weights,
                                   weight_specs$weight_type)
    
    mod_df <- dplyr::bind_rows(mod_df_ls,
                              .id = "dataset_modelling_id")
    
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
  # print(head(mod_df))
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
  
  metr_vec <- c(MuMIn::r.squaredGLMM(m0)[1,], stats::AIC(m0), stats::BIC(m0), MuMIn::AICc(m0))
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
#' This function collects the data needed for the model fitting process, based on 
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
                                 colnames) |>
    unlist() |>
    unique() 
  between_site_resp <- purrr::map(obj$diversity_between$Data,
                                  colnames) |>
    unlist() |>
    unique() 
  
  if (mod_resp %in% within_site_resp){
    
    # If not all covars are in obj$sites, as in a possible phylogenetically corr
    # brms model, this will fail!!
    # Could subset to those covars available in sites - as in INLA
    # And throw a warning?!
    # TODO: - only include available columns!!!
    # 
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

    # mod_df <- obj$sites[,unique(c(mod_covars, c("SS", "SSB", "SSBS")))]
    # find relevant obj$diversity$Data[[i]] 

    # TODO: BUG, Shannon_diversity, and Simpson_diversity also match HillXXX
    # Need to check for exact matches, not grepl matches
    # grepl(paste0("^", resp, "$"), tmp) ??
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
    
  } else if (mod_resp %in% between_site_resp){

    all_covars <- unique(c(mod_covars, c("SSBS1", "SSBS2") ))
    missing_covars <- all_covars[!all_covars %in% colnames(obj$site_contrasts)]
    if (length(missing_covars) > 0){
      message(crayon::yellow(paste0("Warning: covariates specified in the model ", 
                                    "formula are missing; ",
                                    paste(missing_covars, colapse = ", "),
                                    "\n")))
    }
    all_covars_sub <- all_covars[all_covars %in% colnames(obj$site_contrasts)]
    mod_df <- obj$site_contrasts[,all_covars_sub]

    # mod_df <- obj$site_contrasts[,unique(c(mod_covars, c("SSBS1", "SSBS2", "SS", "SSB")))]
    
    # find relevant obj$diversity$Data[[]] - which one contains the response...
    div_data_idx <- which(unlist(lapply(obj$diversity_between$Data,
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
      dplyr::left_join(obj$diversity_between$Data[[div_data_idx]][,c(mod_resp, "SSBS1", "SSBS2")],
                      by = c("SSBS1", "SSBS2"))
    
  } else {
    message(crayon::magenta("Error: response variable specified is not available.\n"))
    stop()
  }
  
  mod_df <- stats::na.omit(mod_df)
  
  # TODO: FOR SAFETY (esp brms), should convert any chr columns to fct (mod_covars)
  # Do this in prep_model_data
  # chr_cols_idx <- purrr::map(mod_covars, 
  #                        function(x, df = mod_df){
  #                          is.character(df[[x]])}) |> 
  #   unlist()
  # chr_cols <- mod_covars[chr_cols_idx]
  # # print(chr_cols)
  # for (c in chr_cols){
  #   mod_df[[c]] <- as.factor(mod_df[[c]])
  # }

  return(mod_df)
  
}




# Prep data -scale options, relevel factor columns

#' Prepare model data for modelling
#' 
#' @description 
#' This function prepares the data for the model fitting based on user supplied
#' flags indicating the scaling of continuous variables and/or the releveling
#' of categroical covariates
#' @details 
#' Expects a dataframe object...
#' @param df A [`dataframe`].
#' @param mod_covars A [`character`] vector indicating the model covariates.
#' @param covars_to_scale A [`character`] vector indicating the response variables to scale.
#' @param scale_func A [`function`] to scale specified continuous variables.
#' @param covars_to_relevel A [`character`] vector containing categorical variables
#' to relevel.
#' @param base_levels A [`character`] vector containing baseline levels
#' @keywords data processing
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
    
    # drop unused factor levels
    df <- droplevels(df)

    # relevel
    for (i in 1:length(covars_to_relevel)){
      if (is.character(df[[covars_to_relevel[i]]])){
        df[[covars_to_relevel[i]]] <- factor(df[[covars_to_relevel[i]]])
      }
      df[[covars_to_relevel[i]]] <- relevel(df[[covars_to_relevel[i]]],
                                             base_levels[i])
    }
  }
  
  return(df)
  
}

# prep_model_data(tmp_db$sites[,c("LUI", "Sampling_effort")],
#                 c("Sampling_effort"),
#                 scale,
#                 c("LUI"), 
#                 c("Cropland_Intense_use"))


# TODO
#' Calculate weights for use in integrated modelling
#' 
#' @description 
#' This function calculates and assigns weights to rows in the modelling data
#' based on the provided arguments
#' @details 
#' Expects a list of dataframes
#' @param df_ls A list of [`dataframes`], each associated with a biodiversity dataset.
#' @param weights A [`numeric`] vector indicating the 'weights' the user wants to 
#' link to each dataset.
#' @param weight_type A [`character`] specifying the form of weighting. Either 
#' "dataset" or "site". 
#' @keywords data processing
#' @examples 
#' \dontrun{
#'   mod_ls <- weight_model_data(df_ls, weights = c(0.5,0.5), weight_type = "dataset") 
#' } 
#' @noRd
weight_model_data <- function(df_ls,
                              model_spec,
                              weights,
                              weight_type){
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
  
  
  # TODO: check weight_type - NULL, dataset, site
  if (is.null(weight_type)){
    message(crayon::yellow("Warning: weight_type missing, setting to 'dataset'."))
    weight_type <- "dataset"
  }
  weight_type <- match.arg(weight_type, c("dataset", "site"), several.ok = F)  
  
  # TODO: check weights - NULL, numeric,.
  if (is.null(weights)){
    weights = rep(1, length(df_ls))
    message(crayon::yellow(paste0("Warning: weights not specified, using equal weight per ",
                                  weight_type, ".")))
  }
  
  
  # TODO: need to know model type... - 
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


# # %%%%%%
# df_ls_1 <- list(df1 = data.frame(x = runif(100)),
#                 df2 = data.frame(x = runif(200)),
#                 df3 = data.frame(x = runif(50)))
# 
# df_ls_2 <- list(df1 = data.frame(x = runif(100)),
#                 df2 = data.frame(x = runif(400)),
#                 df3 = data.frame(x = runif(60)))
# 
# test_1 <- weight_model_data(df_ls_1,
#                             c(0.5,0.5,0.5),
#                             "dataset")
# 
# 
# test_1.1 <- weight_model_data(df_ls_1,
#                             c(0.5,1,0.5),
#                             "site")
# 
# test_2 <- weight_model_data(df_ls_2,
#                             c(1,0.5,0.5),
#                             "dataset")
# 
# 
# test_2.1 <- weight_model_data(df_ls_2,
#                               c(0.5,1,2),
#                               "site")
# 
# lapply(test_1, function(x){sum(x$data_weight)}) |> unlist() |> sum()
# lapply(test_1, function(x){x$data_weight[1]}) |> unlist() 
# 
# lapply(test_1.1, function(x){sum(x$data_weight)}) |> unlist() |> sum()
# lapply(test_1.1, function(x){x$data_weight[1]}) |> unlist() 
# 
# 
# lapply(test_2, function(x){sum(x$data_weight)}) |> unlist() |> sum()
# lapply(test_2, function(x){x$data_weight[1]}) |> unlist() 
# 
# lapply(test_2.1, function(x){sum(x$data_weight)}) |> unlist() |> sum()
# lapply(test_2.1, function(x){x$data_weight[1]}) |> unlist() 
# c(0.000877193, 0.001754386, 0.003508772 ) / max(c(0.000877193, 0.001754386, 0.003508772 ) ) 
# # %%%%%%
