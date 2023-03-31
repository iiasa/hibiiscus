#' Fit an inla-style mixed-effects model
#' 
#' @description 
#' This function fits an inla(bru) model using the supplied database, formula and other
#' function arguments tuning the model fitting process
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @param form A [`formula`] in the style of inla(bru) - including fixed and random effects.
#' @keywords model, fitting
#' @examples 
#' \dontrun{
#'   mod <- hib_fit_inla(obj, ...)
#' } 
#' @export
hib_fit_inla <- function(obj, 
                         # form,
                         # maybe initially separate formula into response (chr), fixed, random (non-spatial)
                         # response,
                         # fixed_covar,
                         response_fixed_form,
                         random_covar,
                         family = "gaussian",
                         fit_settings,
                         incl_spatial_term = T,
                         input_mesh = NULL,
                         mesh_settings = list(mesh_dim = "2d",
                                              max.edge = c(75,100),
                                              cutoff = 1e-12,
                                              offset = NULL),
                         prior_list){
  # Assertions
  assertthat::assert_that(
    inherits(obj, "Database"), # obj is a hib database
    # form can be a formula,
    # is.character(response),
    # assertthat::assert_that(attr(fixed_covar, "class") == "formula"),
    assertthat::assert_that(attr(response_fixed_form, "class") == "formula"),
    assertthat::assert_that(attr(random_covar, "class") == "formula"),
    is.character(family), # family is character
    is.logical(incl_spatial_term),
    inherits(input_mesh,'inla.mesh') || is.null(input_mesh)
  )
  # TODO:
  # Checks on supplied formulas - rhs only for ranefs
  
  # TODO:
  # and maybe rename random_covar to random_form?
  
  # Add progress messages
  message(crayon::green("Collating data for model fitting...\n"))
  
  mod_resp <- all.vars(response_fixed_form[[2]])
  fixed_covars <- all.vars(response_fixed_form[[3]])
  all_covars <- unique(c(fixed_covars,
                         all.vars(random_covar)))
  
  # mod_resp <- all.vars(as.formula(paste0(response, " ~ 1")))
  # 
  # all_covars <- unique(c(all.vars(fixed_covar),
  #                        all.vars(random_covar)))
  
  # NOTE, not all covars are in sites, (Intercept)
  covars_sub <- all_covars[all_covars %in% colnames(obj$sites)]
  
  # Only need lon/lat if incl_spatial_term is TRUE
  if (incl_spatial_term){
    mod_df <- obj$sites[,unique(c(covars_sub, c("SS", "SSB", "SSBS", "Longitude", "Latitude")))]
  } else {
    mod_df <- obj$sites[,unique(c(covars_sub, c("SS", "SSB", "SSBS")))]
  }
  
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
  
  # Add an intercept... 
  if ("Intercept" %in% all_covars){
    mod_df$Intercept <- 1
  }
  
  
  
  # But need to drop "baseline" level... - only if -1 in fixed_covar
  # all.vars(fixed_covar)
  fixed_covar_form <- formula(delete.response(terms(response_fixed_form)))
  
  # Convert categorical fixed-effects to multiple dummy vars...
  X <- model.matrix(fixed_covar_form, data = mod_df)
  
  drop_baseline_level <- function(model_matr, model_df, fixed_covar_form){
    # Loop through fixed eff covars (other than Intercept)
    # Check length of unique cats and compare to N assoc cols in model matrix
    # If the same, drop the first level of the factors
    covars <- all.vars(fixed_covar_form)
    
    for (t in covars){
      if (t != "Intercept"){
        # Find N cols in model matrix
        cats <- unique(as.character(model_df[[t]]))
        cat_terms <- paste0(t, cats)
        n_matr_terms <- sum(cat_terms %in% colnames(model_matr))
        # Calc n levels in model df
        n_cats <- length(cats)
        if (n_matr_terms == n_cats){
          # Drop first term 
          base_idx <- min(which(colnames(model_matr) %in% cat_terms))
          model_matr <- model_matr[,-base_idx]
        }
      }
    }
    return(model_matr)
  }
  
  if (any(grepl("-1", fixed_covars))){
    X <- drop_baseline_level(X, mod_df, fixed_covar_form)
  }
  
  # Make model matrix a df..
  X <- as.data.frame(X)
  # NOTE: inla can't handle spaces of colons in covar names!!
  colnames(X) <- gsub(" ", ".", colnames(X))
  colnames(X) <- gsub(":", "..", colnames(X))
  

  # Add in random effects...
  X <- cbind(X, mod_df[,all.vars(random_covar)])
  
  
  # if response var is not just a column, need to make response var (e.g., log(x))
  # Make y response var - need to make even if it is just a column....
  # A bit hacky!!!
  mod_df <- mod_df %>%
    dplyr::mutate(y = !!response_fixed_form[[2]] )
    # dplyr::mutate(y = !!as.formula(paste0(response, "~ 1"))[[2]] )
  
  # build mesh- if desired...
  if (incl_spatial_term){
    if (is.null(input_mesh)){
      message(crayon::green("Building mesh...\n"))
      # Create a mesh
      # check that mesh_dim, max.edge, cutoff and offset are all in mesh_settings
      # if any not, print out a message telling user that default settings (x..)
      # are being used
      missing_mesh_specs <- c("mesh_dim", "max.edge", "cutoff", "offset")[!(c("mesh_dim", "max.edge", "cutoff", "offset") %in% names(mesh_settings))]
      if ("mesh_dim" %in% missing_mesh_specs){
        message(crayon::yellow("Mesh dimension (2d or 3d) not provided, using default value: 2d \n"))
        mesh_settings$mesh_dim <- "2d"
      }
      if ("max.edge" %in% missing_mesh_specs){
        # TODO: better default - link to location distances...
        message(crayon::yellow("Mesh max.edge not provided, using default values: c(75,100) \n"))
        mesh_settings$max.edge <- c(75,100)
      }
      if ("cutoff" %in% missing_mesh_specs){
        message(crayon::yellow("Mesh cutoff not provided, using default value: 1e-12 \n"))
        mesh_settings$cutoff <- 1e-12
      }
      if ("offset" %in% missing_mesh_specs){
        message(crayon::yellow("Mesh offset not provided, using default values: c(10,20) \n"))
        mesh_settings$offset <- c(10,20)
      }
      
      # TODO: move this to a function - hib_build_inla_mesh <- function(coords, mesh_settings)??
      coords <- as.matrix(mod_df[,c("Longitude", "Latitude")])
      
      
      if (mesh_settings$mesh_dim == "3d"){
        input_mesh <- inla.mesh.create(loc = coords,
                                       globe = T,
                                       refine = list(max.edge = mesh_settings$max.edge))
      }
      
      if (mesh_settings$mesh_dim == "2d"){
        # Try and set some sensible baseline mesh params if not provided
        input_mesh <- inla.mesh.2d(loc = coords,
                                   max.edge = mesh_settings$max.edge,
                                   cutoff = mesh_settings$cutoff,
                                   offset = mesh_settings$offset)
      }
      
      # input_mesh <- inla.mesh.2d(loc = coords,
      #                            max.edge = c(75, 100))
    } 
    # else {
    #   # or use provided mesh
    #   mesh <- input_mesh
    # }
    
    # Set up spde model
    site.spde <- inla.spde2.matern(input_mesh, alpha = 2)
    
    # TODO: add prior specs...
    
    # Hosts.spde = inla.spde2.pcmatern(mesh = Mesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
    # w.Host <- inla.spde.make.index('w', n.spde = site.spde$n.spde) # making the w
    
    # Create weights (A) matrix
    A <- inla.spde.make.A(input_mesh, 
                          loc = coords)
    s.site <- inla.spde.make.index('s', n.spde = site.spde$n.spde)
    
    # collate model fitting data
    # stack data
    stk <- inla.stack(tag = "tr",
                      data = list(y = mod_df$y),
                      A = list(1, A),
                      effects = list(X = X,
                                     # Intercept = rep(1, nrow(mod_df)),
                                     s = s.site))
  } else {
    # Just use model df
    stk <- inla.stack(tag = "tr",
                      data = list(y = mod_df$y),
                      A = list(1),
                      effects = list(X = X))
  }
  
  # build formula 
  resp_chr <- "y ~ "
  
  # Get subset of covars that aren't factors, keep in fixed_covar - conv to char
  fixed_covar_terms_chr <- all.vars(fixed_covar_form)
  factor_fixed_terms_chr <- fixed_covar_terms_chr[lapply(fixed_covar_terms_chr,
                                                 function(x, df){
                                                   is.factor(df[[x]]) | is.character(df[[x]])
                                                 },
                                                 df = mod_df) %>% unlist()]
  # as.character(fixed_covar)[2]
  # gsub each factor_fixed_terms for ""
  # tidy character version of form to ensure no extra "+"
  if (length(factor_fixed_terms_chr)> 0){
    # fixed_chr <- paste(as.character(fixed_covar)[2],
    #                    paste(factor_fixed_terms_chr, collapse = " - "),
    #                    sep = " - ")
    # Above is not correct!!!
    
    terms_to_drop_idx <- which(attr(terms(fixed_covar_form), "term.labels") %in% factor_fixed_terms_chr)
    
    fixed_chr <- as.character(formula(drop.terms(terms(fixed_covar_form), dropx = terms_to_drop_idx)))[2]
    
    # as.character(drop.terms(terms(fixed_covar), dropx = terms_to_drop_idx))
    
    # Get colnames from X, and drop those cols that are already in fixed_covar/random covar terms 
    fact_covars_to_add <- colnames(X)[!colnames(X) %in% all_covars]
    
    if (length(fact_covars_to_add)>0){
      # Paste remaining colnames to add to continuous covars...
      fact_covars_to_add <- paste(fact_covars_to_add, collapse = " + ")
      fixed_chr <- paste(fixed_chr,
                          fact_covars_to_add, 
                          sep = " + ")
    }
  } else {
    fixed_chr <- as.character(fixed_covar_form)[2]
  }
  
  # Paste on random eff specification
  ranef_chr <- as.character(random_covar)[2]
  
  inla_form <- paste0(resp_chr, fixed_chr)
  inla_form <- paste(inla_form, ranef_chr, sep = " + ")
  
  if (incl_spatial_term){
    # add spatial term
    spat_chr <- "f(s, model = site.spde)"
    inla_form <- paste(inla_form, spat_chr, sep = " + ")
  }

  inla_form <- formula(inla_form)
  
  message(crayon::green("Fitting model (this may take a while)...\n"))
  # fit model
  m0 <- inla(formula = inla_form,
                 family = "gaussian",
                 data=inla.stack.data(stk),
                 # control.fixed = 
                 control.compute = list(dic = TRUE, waic = TRUE),
                 control.predictor = list(A = inla.stack.A(stk)),
                 verbose = T)
  # summary(m_inla)
  # return(m_inla)
  # TODO:
  # create inla_mod class
  message(crayon::green("Collating model details ..."))
  
  n_spp <- obj$diversity$N_species_total[[div_data_idx]]

  # instantiate an inla__mod class
  mod_obj <- hib_inla_model$new(
    db_obj = obj,
    mod_df = mod_df,
    inla_obj = m0,
    # TODO:
    user_formula = formula(paste(as.character(response_fixed_form)[2],
                                 as.character(response_fixed_form)[1],
                                 as.character(response_fixed_form)[3],
                                 "+",
                                 as.character(random_covar)[2], sep = " ")),
    inla_fitted_formula = inla_form,
    inla_incl_spatial_term = incl_spatial_term,
    inla_stack = stk,
    inla_mesh = input_mesh,
    inla_mesh_settings = mesh_settings,
    
  # fixef_df = fixef_df,
  #   metr_vec = metr_vec,
    n_spp = n_spp
  # TODO:

  # family = "gaussian",
  # fit_settings,
  # mesh_settings,
  # prior_list
  )
  
  return(mod_obj)
  
  # %%%%%
  # lmer example

  # fixef_df <- as.data.frame(lme4::fixef(m0))
  # colnames(fixef_df) <- "Coefficient_value"
  # fixef_df <- tibble::rownames_to_column(fixef_df, "Coefficient_name")
  # # print(calc_CIs)
  # if (!is.null(calc_CIs)){
  #   message(crayon::green("Getting confidence intervals ..."))
  #   confint_df <- lme4::confint.merMod(m0, parm = "beta_", method = calc_CIs)
  #   confint_df <- as.data.frame(confint_df)
  #   colnames(confint_df) <- c("loCI", "hiCI")
  #   confint_df$Coefficient_name <- rownames(confint_df)
  #   fixef_df <- dplyr::left_join(fixef_df,
  #                                confint_df)
  # }
  # 
  # metr_vec <- c(MuMIn::r.squaredGLMM(m0)[1,], AIC(m0), BIC(m0), MuMIn::AICc(m0))
  # names(metr_vec) <- c("R2m", "R2c", "AIC", "BIC", "AICc")
  
  # %%%%%
}


# %%%%%%%%
# TESTING
# tmp <- log(Total_abundance) ~ -1 + Intercept + pa + man1 + man2 + use1 + 
#   use2 + f(s, model = spde) + f(Binomial, model = "iid") + 
#   f(Sys_realm, model = "iid")
# tmp1 <- Total_abundance ~ x + y + f(s, model = spde)
# 
# raw_resp <- tmp[[2]]
# raw_resp1 <- tmp1[[2]]
# 
# resp_var <- all.vars(tmp[[2]])
# resp_var1 <- all.vars(tmp1[[2]])
# 
# all(as.character(raw_resp) == resp_var)
# all(as.character(raw_resp1) == resp_var1)
# 

# Need to expand mod_df to include factor expansion...
# Using model.matrix
# But only for fixef in covar...

# Random slope in INLA
# f(region, covariate, model ="iid")
# where first var is grouping var, second is covar...
# equiv to (covar | region) in lmer

# %%%%%%%%

