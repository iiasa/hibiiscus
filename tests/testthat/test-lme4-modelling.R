# Testthat checks for hib_fit_lme4.R, cv, projection

test_that("Check lme4 modelling pipeline works as desired.", {
  tmp_db <- create_hib(#path_species = "inst/extdata/dummy_species_rc.rds",
                       path_species = system.file("extdata", "dummy_species_rc.rds", package = "hibiiscus", mustWork = TRUE),
                       version = "custom",
                       custom_version_name = "dummy_data")
  # tmp_db <- load_hib(#"data/dummy_database.rds"
  #                    system.file("extdata", "dummy_database.rds", package = "hibiiscus", mustWork = TRUE))
  
  tmp_db <- hib_calc_div_within(tmp_db,
                                c("Total_abundance", "Species_richness",
                                  "Mean_abundance", "PIE"))
  
  
  tmp_db <- hib_calc_div_between(tmp_db,
                                 c("baselga_bray_balanced", "euclidean", "theta",
                                   "horn"), 
                                 baseline_lui = "Primary_Minimal_use")
  tmp_db <- hib_calc_site_contrasts(tmp_db,
                                    metrics = c("lui_contrast", 
                                                "geogr_dist_hav", "geogr_dist_geod",
                                                "geogr_midpoint"),
                                    baseline_lui = "Primary_Minimal_use")
  
  # within-site model
  # ignore warnings...
  suppressWarnings(mod_w <- hib_fit_lme4(tmp_db,
                                         log(Total_abundance + 1) ~ LUI + (1 | SS),
                                         calc_CIs = "Wald",
                                         relevel_specs = c("LUI" = "Primary_Minimal_use")
                                         ))
  
  expect_true(inherits(mod_w,
                       "HibModelLME4"))
  
  # between-site model
  suppressWarnings(mod_b <- hib_fit_lme4(tmp_db,
                                         Theta ~ LUI_contrast + Geogr_dist_hav + (1 | SS),
                                         calc_CIs = "Wald",
                                         relevel_specs = c("LUI_contrast" = "Primary_Minimal_use_vs_Primary_Minimal_use")))
  expect_true(inherits(mod_b,
                       "HibModelLME4"))
  
  # mod_w$plot()
  # mod_b$plot()
  
  # cv
  suppressWarnings(cv_w <- hib_mod_cv(mod_w,
                                      tmp_db,
                                      c("Random-2", "Random-5"),
                                      backtransform = function(x){exp(x) - 1}))
  expect_true(inherits(cv_w,
                       "HibCV"))
  
  suppressWarnings(cv_b <- hib_mod_cv(mod_b,
                                      tmp_db,
                                      c("Random-2")))
  expect_true(inherits(cv_b,
                       "HibCV"))
  
  # cv_w$plot()
  # cv_b$plot()
  
  
  # projection
  # simulate some lui/contrast data for models
  # 100x100 grid
  proj_df_w <- expand.grid(x = 1:100,
                           y = 1:100) |>
    as.data.frame() |>
    dplyr::mutate(LUI = sample(as.character(unique(mod_w$data$LUI)),
                               dplyr::n(), T))
  
  proj_df_b <- expand.grid(x = 1:100,
                           y = 1:100) |>
    as.data.frame() |>
    dplyr::mutate(LUI_contrast = sample(as.character(unique(mod_b$data$LUI_contrast)),
                                        dplyr::n(), T),
                  Geogr_dist_hav = 0)
  
  # %%%%
  # AOA TESTING
  # devtools::load_all()
  # aoa_w <- hib_calc_aoa(mod_w, cv_w, proj_df_w)
  suppressWarnings(
    proj_w_aoa <- hib_mod_project(mod_w,
                                proj_df_w,
                                backtransform = function(x){exp(x) - 1},
                                calc_aoa = T,
                                cv_obj = cv_w)
  )
  # %%%%
  
  
  # project for df
  expect_true(inherits(proj_w_aoa,
                       "HibProj"))  
  # Check for aoa_res as a list in proj_w_aoa
  expect_true(is.list(proj_w_aoa$aoa))

  proj_b <- hib_mod_project(mod_b,
                            proj_df_b)
  expect_true(inherits(proj_b,
                       "HibProj")) 
  
  # proj_w$plot()
  # proj_b$plot()
  
  # terra based projection
  proj_rast_w <- proj_df_w |>
    dplyr::mutate(LUI_idx = as.numeric(as.factor(LUI))) |>
    dplyr::select(x, y, LUI_idx, ) |>
    terra::rast() |>
    terra::as.factor()
  
  proj_rast_w <- terra::categories(proj_rast_w,
                                   value = proj_df_w |>
                                     dplyr::mutate(LUI_idx = as.numeric(as.factor(LUI))) |>
                                     dplyr::select(LUI_idx, LUI) |>
                                     dplyr::rename(ID = LUI_idx, category = LUI) |>
                                     unique()
  )
  names(proj_rast_w) <- "LUI"
  
  suppressWarnings(
    proj_w_aoa_ <- hib_mod_project(mod_w,
                                 proj_rast_w,
                                 backtransform = function(x){exp(x) - 1},
                                 calc_aoa = T,
                                 cv_obj = cv_w)
  )
    
  expect_true(inherits(proj_w_aoa_,
                       "HibProj"))  
  # Check for aoa_res as a list in proj_w_aoa_
  expect_true(is.list(proj_w_aoa_$aoa))
  
  proj_rast_b <- proj_df_b |>
    dplyr::mutate(LUI_idx = as.numeric(as.factor(LUI_contrast))) |>
    dplyr::select(x, y, LUI_idx) |>
    terra::rast() |>
    terra::as.factor()
  
  proj_rast_b <- terra::categories(proj_rast_b,
                                   value = proj_df_b |>
                                     dplyr::mutate(LUI_idx = as.numeric(as.factor(LUI_contrast))) |>
                                     dplyr::select(LUI_idx, LUI_contrast) |>
                                     dplyr::rename(ID = LUI_idx, category = LUI_contrast) |>
                                     unique()
  )
  names(proj_rast_b) <- "LUI_contrast"
  
  
  proj_rast_b_ <- proj_df_b |>
    dplyr::select(x, y, Geogr_dist_hav) |>
    terra::rast() 
  
  proj_rast_b <- c(proj_rast_b, proj_rast_b_)
  
  proj_b_ <- hib_mod_project(mod_b,
                             proj_rast_b)
  expect_true(inherits(proj_b_,
                       "HibProj"))  
  
  # proj_w_$plot()
  # proj_b_$plot()
  
  # proj_w_$projection
  # proj_b_$projection
})

