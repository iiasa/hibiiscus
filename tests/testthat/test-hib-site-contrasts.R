# Testthat checks for hib_site_contrasts.R

test_that("Check that hib_site_contrasts works as expected.", {
  tmp_db <- hib_create_db(#path_species = "inst/extdata/dummy_species_rc.rds",
                       path_species = system.file("extdata", "dummy_species_rc.rds", package = "hibiiscus", mustWork = TRUE),
                       version = "custom",
                       custom_version_name = "dummy_data")
  # tmp_db <- load_hib(#"data/dummy_database.rds"
  #                    system.file("extdata", "dummy_database.rds", package = "hibiiscus", mustWork = TRUE))
  
  # -Incorrect object input
  expect_error(hib_calc_site_contrasts("hello",
                                       metrics = c("lui_contrast", 
                                                   "geogr_dist_hav", "geogr_dist_geod",
                                                   "geogr_midpoint"),
                                       baseline_lui = "Primary_Minimal_use"),
               regexp = "obj does not inherit from class HibDatabase")
  
  # -Incorrect baseline LUI
  expect_error(hib_calc_site_contrasts(tmp_db,
                                       metrics = c("lui_contrast", 
                                                   "geogr_dist_hav", "geogr_dist_geod",
                                                   "geogr_midpoint"),
                                       baseline_lui = "Primary"),
               regexp = "baseline_lui")
  
  # -Incorrect metric input
  expect_error(hib_calc_site_contrasts(tmp_db,
                                       metrics = c(5, 10),
                                       baseline_lui = "Primary_Minimal_use"),
               regexp = "metrics is not a character vector")
  
  # Correct specification
  tmp_db <- hib_calc_site_contrasts(tmp_db,
                                    metrics = c("lui_contrast", 
                                                "geogr_dist_hav", "geogr_dist_geod",
                                                "geogr_midpoint"),
                                    baseline_lui = "Primary_Minimal_use")
  expect_true(inherits(tmp_db,
                       "HibDatabase"))
  expect_true(is.data.frame(tmp_db$site_contrasts))
  
})
