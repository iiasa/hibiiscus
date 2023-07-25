# Testthat checks for hib_div_metrics.R

# %%%%%%%%%%
# Within-sample
# %%%%%%%%%%

test_that("Check that within-sample diversity metrics work.", {
  tmp_db <- load_hib("data/dummy_database.rds")
  
  # -Won't work if metric is not supported
  expect_error(hib_calc_div_within(tmp_db,
                                   c(5, 10)),
               regexp = "metrics is not a character vector")
  
  # expect_error(hib_calc_div_within(tmp_db,
  #                                  c("not-a-metric")),
  #              regexp = "Error in match")
  
  within_site_metrics <- c("Total_abundance", "Species_richness",
                           "Mean_abundance", "PIE")
  
  tmp_db <- hib_calc_div_within(tmp_db,
                                within_site_metrics)
  
  expect_true(inherits(tmp_db,
                       "Database"))
  
  expect_true(is.data.frame(tmp_db$diversity_within))
  expect_true(nrow(tmp_db$diversity_within) == length(within_site_metrics))
  
  # tmp_db <- hib_calc_div_within(tmp_db,
  #                               c("Total_abundance", "Species_richness",
  #                                 "Geometric_mean_abundance"))

})



# %%%%%%%%%%
# Between-sample
# %%%%%%%%%%
devtools::load_all()

test_that("Check that between-sample diversity metrics work.", {
  tmp_db <- load_hib("data/dummy_database.rds")
  
  # -Won't work if metric is not supported
  expect_error(hib_calc_div_between(tmp_db,
                                   c(5, 10)),
               regexp = "metrics is not a character vector")
  
  
  between_site_metrics <- c("baselga_bray_balanced", "euclidean", "theta",
                            "horn")
  # unique(tmp_db$sites$LUI)
  tmp_db <- hib_calc_div_between(tmp_db,
                                 between_site_metrics, 
                                 baseline_lui = "Primary_Minimal_use")
  
  expect_true(inherits(tmp_db,
                       "Database"))
  
  expect_true(is.data.frame(tmp_db$diversity_between))
  expect_true(nrow(tmp_db$diversity_between) == length(within_site_metrics))
  
  # tmp_db <- hib_calc_div_between(tmp_db,
  #                                c("euclidean", "theta",
  #                                  "hypervolume", "cosine"), 
  #                                baseline_lui = "Primary_Minimal_use")
})
