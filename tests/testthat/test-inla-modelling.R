# Testthat checks for hib_fit_inla.R, cv, projection

test_that("Check inla modelling pipeline works as desired.", {
  tmp_db <- load_hib("data/dummy_database.rds")
  
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
  suppressWarnings(mod_w <- hib_fit_inla(tmp_db,
                                         log(Total_abundance + 1) ~ -1 + LUI,
                                         ~ f(SS, model = "iid")))
  expect_true(inherits(mod_w,
                       "hib_inla_model"))
  
  # between-site model - fails...-inla error
  # suppressWarnings(mod_b <- hib_fit_inla(tmp_db,
  #                                        Theta ~ LUI_contrast + Geogr_dist_hav,
  #                                        ~ f(SS, model = "iid")))
  # expect_true(inherits(mod_b,
  #                      "hib_inla_model"))
  
  # mod_w$plot()
  # mod_b$plot()
  
  # cv
  suppressWarnings(cv_w <- hib_mod_cv(mod_w,
                                      tmp_db,
                                      c("Random-2"),
                                      backtransform = function(x){exp(x) - 1}))
  expect_true(inherits(cv_w,
                       "hib_cv_res"))
  
  # suppressWarnings(cv_b <- hib_mod_cv(mod_b,
  #                                     tmp_db,
  #                                     c("Random-2")))
  # expect_true(inherits(cv_b,
  #                      "hib_cv_res"))
  
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
  
  # proj_df_b <- expand.grid(x = 1:100,
  #                          y = 1:100) |>
  #   as.data.frame() |>
  #   dplyr::mutate(LUI_contrast = sample(as.character(unique(mod_b$data$LUI_contrast)),
  #                                       dplyr::n(), T),
  #                 Geogr_dist_hav = 0)
  
  
  # project for df, ?and terra rast?
  proj_w <- hib_mod_project(mod_w,
                            proj_df_w,
                            backtransform = function(x){exp(x) - 1})
  expect_true(inherits(proj_w,
                       "hib_mod_projection"))  
  
  # proj_b <- hib_mod_project(mod_b,
  #                           proj_df_b)
  # expect_true(inherits(proj_b,
  #                      "hib_mod_projection")) 
  
  # proj_w$plot()
  # proj_b$plot()
  
  
  # %%%%%
  
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
  
  proj_w_ <- hib_mod_project(mod_w,
                             proj_rast_w,
                             backtransform = function(x){exp(x) - 1})
  expect_true(inherits(proj_w_,
                       "hib_mod_projection"))  
  
  
  # proj_rast_b <- proj_df_b |>
  #   dplyr::mutate(LUI_idx = as.numeric(as.factor(LUI_contrast))) |>
  #   dplyr::select(x, y, LUI_idx) |>
  #   terra::rast() |>
  #   terra::as.factor()
  # 
  # proj_rast_b <- terra::categories(proj_rast_b,
  #                                  value = proj_df_b |>
  #                                    dplyr::mutate(LUI_idx = as.numeric(as.factor(LUI_contrast))) |>
  #                                    dplyr::select(LUI_idx, LUI_contrast) |>
  #                                    dplyr::rename(ID = LUI_idx, category = LUI_contrast) |>
  #                                    unique()
  # )
  # names(proj_rast_b) <- "LUI_contrast"
  # 
  # 
  # proj_rast_b_ <- proj_df_b |>
  #   dplyr::select(x, y, Geogr_dist_hav) |>
  #   terra::rast() 
  # 
  # proj_rast_b <- c(proj_rast_b, proj_rast_b_)
  # 
  # proj_b_ <- hib_mod_project(mod_b,
  #                            proj_rast_b)
  # expect_true(inherits(proj_b_,
  #                      "hib_mod_projection"))  
})
