# Testthat checks for update_hib.R

test_that("Check that update_hib works as desired.", {
  
  tmp_db_1 <- load_hib("data/dummy_database.rds")
  
  # Create a second db from random simul func - saved...
  # tmp_spp_2 <- create_dummy_data("species", 15)
  # saveRDS(tmp_spp_2, "data/dummy_species_rc_2.rds")
  
  # -Doesn't work if ....
  # -Path to new data is missing
  expect_error(update_hib(tmp_db_1, 
                          path_species = NULL,
                          version = "custom",
                          custom_version_name = "dummy_data_2"),
               regexp = "path_species is not a character vector")
  # -Path to new data is wrong
  expect_error(update_hib(tmp_db_1, 
                          path_species = "data/not_real_dummy_data.rds",
                          version = "custom",
                          custom_version_name = "dummy_data_2"),
               regexp = "Path 'data/not_real_dummy_data.rds' does not exist")
  # -Original object is not a hib database
  expect_error(update_hib("not a database", 
                          path_species = "data/not_real_dummy_data.rds",
                          version = "custom",
                          custom_version_name = "dummy_data_2"),
               regexp = "obj does not inherit from class Database")
  
  
  # Does work
  tmp_db_upd <- update_hib(tmp_db_1, 
                           path_species = "data/dummy_species_rc_2.rds",
                           version = "custom",
                           custom_version_name = "dummy_data_2")
  expect_true(inherits(tmp_db_upd,
                       "Database"))
  # tmp_db_upd$version
  # 
})
