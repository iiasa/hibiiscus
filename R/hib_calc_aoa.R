#' Calculate an Area-of-Applicability (AOA) given a model, cross-validation output, and
#' the data being used to make projections
#' 
#' @description
#' This function estimates an AOA, for each of the cross-validation approaches, 
#' following the methods of Meyer and Pebesma 2021
#' @details
#' Expects a [`hib_X_model`] object, a [`hib_cv_res`] object, and data used to 
#' make a projection.
#' @param mod_obj A [`hib_X_model`] object.
#' @param cv_obj A [`hib_cv_res`] object, based on the above mod_obj.
#' @param proj_data A [`data.frame`], [`SpatRaster`] or list of [`SpatRaster`], 
#' containing covariates for projecting the model.
#' @keywords model, projections, 
#' @examples 
#' \dontrun{
#'   aoa_res <- hib_calc_aoa(mod_obj, cv_obj, proj_data)
#' } 
#' @export
hib_calc_aoa <- function(mod_obj, 
                         cv_obj,
                         proj_data){
  # Note: proj_data should not include prediction - just covars...
  
  assertthat::assert_that(
    inherits(mod_obj, "hib_lme4_model") || inherits(mod_obj, "hib_inla_model") || inherits(mod_obj, "hib_brms_model"), # obj is a hib model
    inherits(cv_obj, "hib_cv_res"),
    is.data.frame(proj_data) || inherits(proj_data, "SpatRaster") || ( is.list(proj_data) & !is.data.frame(proj_data) & all(unlist(lapply(proj_data, function(x){inherits(x, "SpatRaster")}))) )
  )
  

  # Get covars and response names
  # lme4:
  if (inherits(mod_obj, "hib_lme4_model")){
    all_covars <- all.vars(formula(delete.response(terms(formula(mod_obj$model)))))
    rand_vars <- names(lme4::ranef(mod_obj$model))
    fix_vars <- all_covars[!all_covars %in% rand_vars]
    
    resp_form <- formula(mod_obj$model)[[2]]
  }
  
  # inla:
  if (inherits(mod_obj, "hib_inla_model")){
    # Package check
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop(
        "Package \"INLA\" must be installed to use this function.",
        call. = FALSE
      )
    }
    # OR
    # rlang::check_installed("INLA", reason = "to use INLA-based model fitting.")
    all_covars <- all.vars(formula(delete.response(terms(mod_obj$user_formula$full))))
    fix_vars <- all.vars(formula(delete.response(terms(mod_obj$user_formula$response_fixed))))
    # ignore "Intercept" too
    fix_vars <- fix_vars[!fix_vars == "Intercept"]
    
    resp_form <- mod_obj$user_formula$response_fixed[[2]]

    # print(fix_vars)
    # print(resp_form)
  }
  
  # brms:
  if (inherits(mod_obj, "hib_brms_model")){
    all_covars <- all.vars(formula(delete.response(terms(mod_obj$model$formula$formula))))
    rand_vars <- names(brms::ranef(mod_obj$model))
    fix_vars <- all_covars[!all_covars %in% rand_vars]
    # Intercept??
    fix_vars <- fix_vars[!fix_vars == "Intercept"]
    
    resp_form <- mod_obj$model$formula$formula[[2]]
  }
  
  
  # Proceed with DI calculation process
  
  # TRAIN DI
  
  # 1. fit random forest to all data - need to use transformed response...
  mod_df <- mod_obj$data
  mod_df <- mod_df |>
    dplyr::mutate(y = !!resp_form )
  
  mod_df <- mod_df[, c(fix_vars, "y"), drop = F]
  
  m_rf <- caret::train(mod_df[, fix_vars, drop = F ],
                       mod_df$y,
                       method="rf",
                       importance=TRUE,
                       trControl = caret::trainControl(method="cv"))
  
  # 2. get weights... - 
  
  # weight <- CAST:::aoa_get_weights(m_rf, variables = fix_vars)
  weight <- aoa_get_weights_local(m_rf, fix_vars)
  
  # 3. convert cat data to dummy (can use aoa_categorial_train(train, variables, weight = 1)
  # use own approach)
  # catupdate <- CAST:::aoa_categorial_train(mod_df[, fix_vars, drop = F ],
  #                                          variables = fix_vars, 
  #                                          weight = weight)
  # aoa_categorial_train_local
  catupdate <- aoa_categorial_train_local(mod_df[, fix_vars, drop = F ],
                                          variables = fix_vars, 
                                          weight = weight)
  

  train <- catupdate$train
  weight <- catupdate$weight
  
  
  # 4. scale data
  # print(head(train))
  train <- scale(train[, names(weight)])
  
  # make sure all variables have variance
  # if (any(apply(train, 2, FUN=function(x){all(is.na(x))}))){
  #   stop("some variables in train seem to have no variance")
  # }
  # 
  # save scale param for later
  scaleparam <- attributes(train)
  
  # 4.1 multiply covar by weights
  # make weight an array, dims equiv to train, byrow = T
  weight_matr <- matrix(unlist(weight), 
                        ncol = ncol(train),
                        nrow = nrow(train),
                        byrow = T)
  
  train_mult <- train[, names(weight)] * weight_matr
  
  train_mult <- as.data.frame(train_mult)
  colnames(train_mult) <- names(weight)
  
  train_mult$SSBS <- mod_obj$data$SSBS
  
  # return(train_mult)
  
  # 5. DI for each cv-block for each cv-type - need to link to SSBS/fold...
  N_cv_type <- length(cv_obj$cv_results)
  
  DI_ls <- vector("list", N_cv_type)
  
  
  for (i in 1:N_cv_type){
    # for each cv_type
    # get fold_id
    tmp <- train_mult |>
      dplyr::left_join(cv_obj$cv_results[[i]]$cv_pred_df[, c("SSBS", "fold_id")],
                       by = "SSBS")
    
    # fold_type:
    tmp_fold_type <- cv_obj$cv_results[[i]]$fold_type
    
    # for each fold, 
    tmp_fold_ids <- unique(tmp$fold_id)

    N_folds <- NULL
    N_folds <- length(tmp_fold_ids)
    tmp_di <- vector("list", N_folds)
    
    for (j in 1:N_folds){
      # subset to ssbs of interest, calc DI - euclidean
      tmp_sub <- tmp |>
        dplyr::filter(fold_id == tmp_fold_ids[j])
      
      # get min/mean of rows
      tmp_dist <- tmp_sub |> 
        dplyr::select(-c(SSBS, fold_id)) |>
        dist() |>
        as.matrix()
      
      diag(tmp_dist) <- NA
      
      tmp_sub$row_mean <- rowMeans(tmp_dist, na.rm = T)
      tmp_sub$row_min <- matrixStats::rowMins(tmp_dist, na.rm = T)
      
      tmp_di[[j]] <- tmp_sub[, c("SSBS", "fold_id", "row_min", "row_mean")]
    }
    # normalise DI by mean of row means
    tmp_di_df <- dplyr::bind_rows(tmp_di)
    # 6. normalise DI by avg
    tmp_di_df <- tmp_di_df |>
      dplyr::mutate(DI = row_min / mean(row_mean, na.rm = T))
    
    # 7. threshold DI for each cv-type
    DI_thr <- quantile(tmp_di_df$DI, 0.75, na.rm = T) + (1.5 * stats::IQR(tmp_di_df$DI, na.rm=T))
    
    # add info to list (cv_type, DI_thr, DI_mean_row_mean, DI_df, )
    DI_ls[[i]] <- list(fold_type = tmp_fold_type,
                       DI_thr    = DI_thr,
                       DI_mean_row_mean = mean(tmp_di_df$row_mean, na.rm = T),
                       DI_df     = tmp_di_df)
    
  }
  
  
  # TEST DI
  # 8. return DI for proj data..., indicate aoa at different cv-types - 
  #    need to scale proj based on train scale... calc dist between proj and train...
  # Could just use mod_df??
  mod_data <- mod_obj$data
  # Need to think about proj format, if df, easy - as above...
  if (is.data.frame(proj_data) || inherits(proj_data, "SpatRaster")){
    proj_aoa <- calc_test_DI(DI_ls, 
                             train_mult,
                             weight,
                             scaleparam,
                             mod_data,
                             fix_vars,
                             proj_data)
  }
  if (( is.list(proj_data) & 
        !is.data.frame(proj_data) & 
        all(unlist(lapply(proj_data, function(x){inherits(x, "SpatRaster")}))) )){
    proj_aoa <- purrr::map(proj_data,
                           function(p_dat,
                                    tr_DI = DI_ls,
                                    tr_mult = train_mult,
                                    var_w = weight,
                                    s_param = scaleparam,
                                    m_dat = mod_data,
                                    fix_v = fix_vars){
                             calc_test_DI(tr_DI, 
                                          tr_mult,
                                          var_w,
                                          s_param,
                                          m_dat,
                                          fix_v,
                                          p_dat)})
  }
  
  return(list(variable_weights = weight,
              variable_scale_params = scaleparam,
              # To output - expanded, unscaled, un weighted train and proj??
              # train_dat = train,
              train_DI = DI_ls,
              # proj_dat = proj_aoa$proj_mult,
              proj_aoa = proj_aoa$proj_aoa))
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal helper functions

 
#' Calculate a Dissimilarity Index (DI) between the data used to fit the model,
#' and the data being used for the projection
#' 
#' @description
#' This function estimates DI for each data point in the projection data and under 
#' each of the cross-validation approaches.
#' -following the methods of Meyer and Pebesma 2021
#' @details
#' Expects information on DI for each cv approach, model-specific data, and 
#' data used to make a projection.
#' @param train_DI A [`list`] containing information on the DI values for "training"
#' data under a range of cross-validation strategies
#' @param train_mult A [`data.frame`] containing the model fitting data, weighted 
#' by variable importance
#' @param variable_weights A [`data.frame`] containing variable importance weights
#' @param variable_scale_params A [`list`] of the parameters used to scale each variable
#' @param mod_data A [`data.frame`] of the data used to fit the original model
#' @param fix_vars A [`character`] identifying model fixed effects
#' @param proj_data A [`dataframe`], [`SpatRaster`] or list of [`SpatRaster`], 
#' containing covariates for projecting the model.
#' @keywords projections, dissimilarity
#' @examples 
#' \dontrun{
#'   di_out <- calc_test_DI(train_DI, train_mult, variable_weights,
#'                          variable_scale_params, mod_data, fix_vars, proj_data)
#' } 
#' @noRd
calc_test_DI <- function(train_DI, 
                         train_mult,
                         variable_weights,
                         variable_scale_params,
                         mod_data,
                         fix_vars,
                         proj_data){
  assertthat::assert_that(
    is.list(train_DI),
    is.data.frame(train_mult) || is.matrix(train_mult),
    is.data.frame(variable_weights),
    is.list(variable_scale_params),
    is.data.frame(mod_data),
    is.character(fix_vars),
    is.data.frame(proj_data) || inherits(proj_data, "SpatRaster")
  )
  
  if (is.data.frame(proj_data)){
    data_df <- proj_data
  } 
  if (inherits(proj_data, "SpatRaster")){
    # Convert factor cols to numeric...
    if (any(is.factor(proj_data))){
      proj_data[[which(is.factor(proj_data))]] <- as.numeric(proj_data[[which(is.factor(proj_data))]])
    }
    
    data_df <- terra::as.data.frame(proj_data, xy = T)
  }
  
  # If fix vars are all in data_df, expand categorical covars
  # cols_in_fix_vars <- colnames(data_df)[colnames(data_df) %in% fix_vars]
  # cols_not_in_fix_vars <- colnames(data_df)[!colnames(data_df) %in% fix_vars]
  # 
  fix_vars_in_df <- fix_vars[fix_vars %in% colnames(data_df)]
  # For each fix/cat_var, expand to levels using caret - 
  cat_vars <- fix_vars_in_df[which(sapply(data_df[, fix_vars_in_df, drop = F], class) %in% c("factor","character"))]
  
  # print(cat_vars)
  
  if (length(cat_vars) > 0){
    for (cat_var in cat_vars){
      # print(cat_var)
      data_df[, cat_var] <- factor(as.character(data_df[, cat_var, drop = T]),
                                   levels = levels(mod_data[, cat_var, drop = T]))
      
      data_df[!data_df[, cat_var] %in% unique(mod_data[, cat_var, drop = T]), cat_var] <- NA
      # 
      data_df[, cat_var] <- droplevels(data_df[, cat_var])
      
      # print(unique(mod_data[, cat_var]))
      # print(unique(data_df[, cat_var]))
      
      # then create dummy variables for the remaining levels in train:
      data_dummy <- predict(caret::dummyVars(paste0("~", cat_var), data = mod_data, 
                                             sep = ""), data_df)
      data_dummy[is.na(data_df[, cat_var]),] <- 0
      
      data_df <- data.frame(data_df, data_dummy)
    }
    data_df <- data_df[, -which(names(data_df) %in% cat_vars)]
  }
  
  # print("HERE")
  # For each non-fix var, get those that match names in weights...
  # Having expanded /collated covar-level names, check all are present (compare to names(weight))
  
  missing_vars <-  names(variable_weights)[!names(variable_weights) %in% colnames(data_df)]
  if (length(missing_vars) > 0){
    message(crayon::red(paste0("Error: variables missing from projection data; ",
                               paste(missing_vars, collapse = ", ") ) ))
    stop("Missing covariate, see above.", call. = F)
  }
  
  assertthat::assert_that(
    all(names(variable_weights) %in% colnames(data_df)) 
  )
  # select those cols for scaling - names(weight) 
  # - ensure scaleparams and columsn are in correct order
  data_df[, names(variable_weights)] <- scale(data_df[, names(variable_weights)],
                                    center = variable_scale_params$`scaled:center`[names(variable_weights)],
                                    scale  = variable_scale_params$`scaled:scale`[names(variable_weights)])
  
  # print(variable_scale_params$`scaled:scale`)
  # print(head(data_df))
  
  weight_matr_ <- matrix(unlist(variable_weights), 
                         ncol = ncol(data_df[, names(variable_weights)]),
                         nrow = nrow(data_df),
                         byrow = T)
  
  data_mult <- data_df[, names(variable_weights)] * weight_matr_
  
  data_mult <- as.data.frame(data_mult)
  colnames(data_mult) <- names(data_mult)
  
  data_mult$x <- data_df$x
  data_mult$y <- data_df$y
  
  data_dist_df <- data_mult[, c("x", "y")]

  data_dist_df$min_dist_raw <- c(FNN::knnx.dist(train_mult[, names(variable_weights)], 
                                                data_mult[, names(variable_weights)], 
                                                1))
  
  
  # Calc DI 
  # For each train_DI list item, 
  for (x in train_DI){
    #   add cols to data_dist_df 
    data_dist_df[, paste0("DI_", x$fold_type)] <- data_dist_df$min_dist_raw / x$DI_mean_row_mean
    data_dist_df[, paste0("AOA_", x$fold_type)] <- as.numeric(data_dist_df[, paste0("DI_", x$fold_type)] <= x$DI_thr)
  }
  
  data_dist_df <- data_dist_df |>
    dplyr::select(-min_dist_raw)

  # convert data_dist_df back to proj_data type... - i.e. SpatRast
  if (inherits(proj_data, "SpatRaster")){
    out <- terra::rast(data_dist_df,
                       crs = terra::crs(proj_data))
    # Check dimensions of proj_dat v out
    assertthat::assert_that(
      terra::nrow(out) == terra::nrow(proj_data),
      terra::ncol(out) == terra::ncol(proj_data),
      terra::crs(out) == terra::crs(proj_data)
    )
  } else {
    out <- data_dist_df
    # Check dimensions of proj_dat v out
    assertthat::assert_that(
      nrow(proj_data) == nrow(out)
    )    
  }

  # Also return scaled, multiplied covariate values?? - coult be useful

  return(list(proj_mult = data_mult,
              proj_aoa  = out))
}


#' Collect variable importance weights from a fitted caret train object
#' 
#' @description
#' This function is a local version of the function from the CAST package
#' -following the methods of Meyer and Pebesma 2021
#' - key difference is inclusion of ", drop = F" in weight collection
#' @keywords projections, dissimilarity, weights
#' @noRd
aoa_get_weights_local <- function(model, variables){
  
  assertthat::assert_that(
    inherits(model, "train"),
    is.character(variables)
  )
  weight <- tryCatch(if(model$modelType=="Classification"){
    as.data.frame(t(apply(caret::varImp(model,scale=F)$importance,1,mean)))
  }else{
    as.data.frame(t(caret::varImp(model,scale=F)$importance[,"Overall"]))
  }, error=function(e) e)
  if(!inherits(weight, "error") & length(variables)>1){
    names(weight)<- rownames(caret::varImp(model,scale=F)$importance)
  }else{
    # set all weights to 1
    weight <- as.data.frame(t(rep(1, length(variables))))
    names(weight) = variables
    message("note: variables were not weighted either because no weights or model were given,
    no variable importance could be retrieved from the given model, or the model has a single feature.
    Check caret::varImp(model)")
  }
  
  #set negative weights to 0
  if(!inherits(weight, "error")){
    weight <- weight[, stats::na.omit(match(variables, names(weight))), drop = F]
    if (any(weight<0)){
      weight[weight<0]<-0
      message("negative weights were set to 0")
    }
  }
  return(weight)
  
}


#' Convert categorical/factor variables to expanded, dummy variables for subsequent
#' DI calculations
#' 
#' @description
#' This function is a local version of the function from the CAST package
#' - following the methods of Meyer and Pebesma 2021
#' - key difference is inclusion of 'sep = ""' in variable expansion
#' @keywords projections, dissimilarity, covariates
#' @noRd
aoa_categorial_train_local <- function(train, variables, weight){
  
  assertthat::assert_that(
    is.data.frame(train),
    is.character(variables),
    is.list(weight)
  )
  # get all categorial variables
  catvars <- tryCatch(names(train)[which(sapply(train[,variables, drop = F], class)%in%c("factor","character"))],
                      error=function(e) e)
  
  if (!inherits(catvars,"error")&length(catvars)>0){
    message("warning: predictors contain categorical variables. The integration is currently still under development. Please check results carefully!")
    
    for (catvar in catvars){
      # print(catvar)
      # mask all unknown levels in newdata as NA (even technically no predictions can be made)
      train[,catvar]<-droplevels(train[,catvar])
      
      # print(head(train[,catvar]))

      # then create dummy variables for the remaining levels in train:
      dvi_train <- predict(caret::dummyVars(paste0("~",catvar), data = train,
                                            sep = ""),
                           train)
      train <- data.frame(train,dvi_train)
      # print(head(dvi_train))
      if(!inherits(weight, "error")){
        addweights <- data.frame(t(rep(weight[,which(names(weight)==catvar)],
                                       ncol(dvi_train))))
        names(addweights)<- colnames(dvi_train)
        weight <- data.frame(weight,addweights)
      }
    }
    if(!inherits(weight, "error")){
      weight <- weight[,-which(names(weight)%in%catvars)]
    }
    train <- train[,-which(names(train)%in%catvars)]
  }
  return(list(train = train, weight = weight, catvars = catvars))
}
