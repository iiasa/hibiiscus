#' Fit an inla-style mixed-effects model
#' 
#' @description 
#' This function fits an inla(bru) model using the supplied database, formula and other
#' function arguments tuning the model fitting process
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @param response_fixed_form A [`formula`] with the response and fixed effects
#' @param random_covar  A [`formula`] with random, non-spatial, iid effects specified
#' @param family A [`character`] indicating the model family, e.g., "gaussian"
#' @param fit_settings In development
#' @param incl_spatial_term A [`logical`], whether or not to include a spatial 
#' random effect
#' @param input_mesh Optionally provide an inla mesh object
#' @param mesh_settings Optionally specify inla mash settings, e.g., max.edge
#' @param prior_list In development
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
    assertthat::assert_that(attr(random_covar, "class") == "formula"), # could be null??
    is.character(family), # family is character
    is.logical(incl_spatial_term),
    inherits(input_mesh,'inla.mesh') || is.null(input_mesh),
    is.list(mesh_settings)
  )
  # TODO:
  # maybe rename random_covar to random_form?
  
  # Add progress messages
  message(crayon::green("Collating data for model fitting...\n"))
  
  mod_resp <- all.vars(response_fixed_form[[2]])
  fixed_covars <- all.vars(response_fixed_form[[3]])
  # Could fix_vars be empty...- need at least global intercept...
  all_covars <- unique(c(fixed_covars,
                         all.vars(random_covar)))
  
  
  # TODO
  # have single form input, resp ~ fixed + random   --- ignoring spatial
  
  
  # mod_resp <- all.vars(as.formula(paste0(response, " ~ 1")))
  # 
  # all_covars <- unique(c(all.vars(fixed_covar),
  #                        all.vars(random_covar)))
  
  
  # %%%%%
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

    # NOTE, not all covars are in sites, (Intercept)
    covars_sub <- all_covars[all_covars %in% colnames(obj$sites)]
    
    # Only need lon/lat if incl_spatial_term is TRUE
    if (incl_spatial_term){
      mod_df <- obj$sites[,unique(c(covars_sub, c("SS", "SSB", "SSBS", "Longitude", "Latitude")))]
    } else {
      mod_df <- obj$sites[,unique(c(covars_sub, c("SS", "SSB", "SSBS")))]
    }
    
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
    
    # NOTE, not all covars are in sites, (Intercept)
    covars_sub <- all_covars[all_covars %in% colnames(obj$site_contrasts)]
    
    # Only need lon/lat if incl_spatial_term is TRUE
    if (incl_spatial_term){
      if (("Latitude_midpoint" %in% colnames(obj$site_contrasts)) &
          ("Longitude_midpoint" %in% colnames(obj$site_contrasts))){
        mod_df <- obj$site_contrasts[,unique(c(covars_sub, c("SSBS1", "SSBS2", "SS", "SSB", 
                                                             "Longitude_midpoint", 
                                                             "Latitude_midpoint")))]
      } else {
        message(crayon::yellow("Warning: coordinate inforamtion unavailable, spatial random effect cannot be fitted.\n"))
        mod_df <- obj$site_contrasts[,unique(c(covars_sub, c("SSBS1", "SSBS2", "SS", "SSB")))]
        incl_spatial_term <- F
      }

    } else {
      mod_df <- obj$site_contrasts[,unique(c(covars_sub, c("SSBS1", "SSBS2", "SS", "SSB")))]
    }
    
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
    
    # print(head(mod_df))
    
    n_spp <- obj$diversity_between$N_species_total[[div_data_idx]]
    
  } else {
    message(crayon::magenta("Error: response variable specified is not available.\n"))
    stop()
  }
  # %%%%%

  
  # Add an intercept... 
  if ("Intercept" %in% all_covars){
    mod_df$Intercept <- 1
  }
  
  
  # if response var is not just a column, need to make response var (e.g., log(x))
  # Make y response var - need to make even if it is just a column....
  # A bit hacky!!!
  mod_df <- mod_df %>%
    dplyr::mutate(y = !!response_fixed_form[[2]] )
    # dplyr::mutate(y = !!as.formula(paste0(response, "~ 1"))[[2]] )
  
  X <- mod_df %>%
    dplyr::select(all_covars)
  
  # build mesh- if desired...
  if (incl_spatial_term){
    if (is.null(input_mesh)){
      message(crayon::green("Building mesh...\n"))
      
      if (("Longitude" %in% colnames(mod_df)) &
          ("Latitude" %in% colnames(mod_df)) ){
        coords <- as.matrix(mod_df[,c("Longitude", "Latitude")])
      }
      if (("Longitude_midpoint" %in% colnames(mod_df)) &
          ("Latitude_midpoint" %in% colnames(mod_df)) ){
        coords <- as.matrix(mod_df[,c("Longitude_midpoint", "Latitude_midpoint")])
      }
      
      # Create a mesh
      # check that mesh_dim, max.edge, cutoff and offset are all in mesh_settings
      # if any not, print out a message telling user that default settings (x..)
      # are being used
      missing_mesh_specs <- c("mesh_dim", "max.edge", "cutoff", "offset")[!(c("mesh_dim", "max.edge", "cutoff", "offset") %in% names(mesh_settings))]
      if ("mesh_dim" %in% missing_mesh_specs){
        message(crayon::yellow("Mesh dimension (2d or 3d) not provided, using default value: 2d \n"))
        mesh_settings$mesh_dim <- "2d"
      }
      
      # TODO: better defaults - link to location distances...
      if ( ("max.edge" %in% missing_mesh_specs) &
          (mesh_settings$mesh_dim == "2d") ){
        message(crayon::yellow("Mesh max.edge not provided, using default values: c(75,100) \n"))
        mesh_settings$max.edge <- c(75,100)
        # From ibis:
        # A good guess here is usally a max.edge of between 1/3 to 1/5 of the spatial range.
        # c(diff(range(coords[,1], na.rm = T))/(3*5) , diff(range(coords[,1], na.rm = T))/(3*5) * 2)
      }
      if ( ("cutoff" %in% missing_mesh_specs) &
           (mesh_settings$mesh_dim == "2d") ){
        message(crayon::yellow("Mesh cutoff not provided, using default value: 1e-12 \n"))
        mesh_settings$cutoff <- 1e-12
        # From ibis:
        # Specify as minimum distance between y coordinates
        # val <- min(abs(diff(coords[,2])), na.rm = T)
        # cutoff <- ifelse(val == 0, 1e-12, val)
      }
      if ( ("offset" %in% missing_mesh_specs) &
           (mesh_settings$mesh_dim == "2d") ){
        message(crayon::yellow("Mesh offset not provided, using default values: c(10,20) \n"))
        mesh_settings$offset <- c(10,20)
        # From ibis:
        # Specify offset as 1/100 of the boundary distance
      }
      
      # TODO: move this to a function - hib_build_inla_mesh <- function(coords, mesh_settings)??
      
      if (mesh_settings$mesh_dim == "3d"){
        coords <- INLA::inla.mesh.map(coords, projection = "longlat")
        input_mesh <- INLA::inla.mesh.create(loc = coords,
                                       globe = 8)
        # globe = T,
        # refine = list(max.edge = mesh_settings$max.edge)
      }
      
      if (mesh_settings$mesh_dim == "2d"){
        # Try and set some sensible baseline mesh params if not provided
        input_mesh <- INLA::inla.mesh.2d(loc = coords,
                                   max.edge = mesh_settings$max.edge,
                                   cutoff = mesh_settings$cutoff,
                                   offset = mesh_settings$offset)
      }
      
    } 

    
    # Set up spde model
    site.spde <- INLA::inla.spde2.matern(input_mesh, alpha = 2)
    
    # TODO: add prior specs...
    # TODO: pcmatern ??
    
    # Hosts.spde = inla.spde2.pcmatern(mesh = Mesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
    # w.Host <- inla.spde.make.index('w', n.spde = site.spde$n.spde) # making the w
    
    # Create weights (A) matrix
    A <- INLA::inla.spde.make.A(input_mesh, 
                          loc = coords)
    s.site <- INLA::inla.spde.make.index('s', n.spde = site.spde$n.spde)
    
    # collate model fitting data
    # stack data
    stk <- INLA::inla.stack(tag = "fit",
                      data = list(y = mod_df$y),
                      A = list(1, A),
                      effects = list(X = X,
                                     # Intercept = rep(1, nrow(mod_df)),
                                     s = s.site))
  } else {
    # Just use model df
    stk <- INLA::inla.stack(tag = "fit",
                      data = list(y = mod_df$y),
                      A = list(1),
                      effects = list(X = X))
  }
  
  # build formula 
  resp_chr <- "y ~ "
  
  fixed_covar_form <- formula(delete.response(terms(response_fixed_form)))
  fixed_chr <- as.character(fixed_covar_form)[2]
  
  inla_form <- paste0(resp_chr, fixed_chr)
  
  # Paste on random eff specification
  ranef_chr <- as.character(random_covar)[2]
  inla_form <- paste(inla_form, ranef_chr, sep = " + ")
  
  if (incl_spatial_term){
    # add spatial term
    spat_chr <- "f(s, model = site.spde)"
    inla_form <- paste(inla_form, spat_chr, sep = " + ")
  }

  
  inla_form <- formula(inla_form)

  message(crayon::green("Fitting model (this may take a while)...\n"))
  
  # TODO: if fixed effects, incl a while to ensure fixed effects are fitted...
  # %%%%%
  # fit model
  
  if (length(fixed_covars) > 0){
    check_flag <- 0
    while(check_flag == 0){
      m0 <- INLA::inla(formula = inla_form,
                 family = "gaussian",
                 data = INLA::inla.stack.data(stk),
                 control.fixed = list(expand.factor.strategy='inla'),
                 # control.inla = list(strategy = "", int.strategy = ""), 
                 control.compute = list(dic = TRUE, waic = TRUE),
                 control.predictor = list(A = INLA::inla.stack.A(stk)),
                 verbose = F)

      check_flag <- nrow(m0$summary.fixed)
    }
  } else {
    m0 <- INLA::inla(formula = inla_form,
               family = "gaussian",
               data = INLA::inla.stack.data(stk),
               control.fixed = list(expand.factor.strategy='inla'),
               # control.inla = list(strategy = "", int.strategy = ""), 
               control.compute = list(dic = TRUE, waic = TRUE),
               control.predictor = list(A = INLA::inla.stack.A(stk)),
               verbose = F)
  }

  message(crayon::green("Collating model details ..."))
  
  # instantiate an inla_mod class
  mod_obj <- hib_inla_model$new(
    db_obj = obj,
    mod_df = mod_df,
    inla_obj = m0,
    # TODO:
    user_formula = list(full = formula(paste(as.character(response_fixed_form)[2],
                                 as.character(response_fixed_form)[1],
                                 as.character(response_fixed_form)[3],
                                 "+",
                                 as.character(random_covar)[2], sep = " ")),
                        response_fixed = response_fixed_form,
                        random = random_covar
                        ),
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


