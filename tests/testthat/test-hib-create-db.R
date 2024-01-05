# Tests for hib_create_db.R

# %%%%%%%%%%
# create_hib 
# %%%%%%%%%%


test_that("Check that hib_create_db correclty loads data.", {
  # - check that passing an incorrect file path returns an error
  expect_error(hib_create_db("inst/extdata/not_a_file.rds"),
               regexp = "does not exist")
  # - numeric path returns an error
  expect_error(hib_create_db(45),
               regexp = "not a character vector")
  # # - ensure that incorrect dataset_name returns error
  # expect_error(hib_create_db(#"inst/extdata/dummy_species_rc.rds",
  #                         path_species = system.file("extdata", "dummy_species_rc.rds", package = "hibiiscus", mustWork = TRUE),
  #                         dataset_name = "dummy_data"),
  #              regexp = "'arg' should be one of \"2016\", \"2022\", \"custom\"")
  
  # - check that only if correct files, and dataset_name are included will it work 
  tmp_db <- hib_create_db(#path_species = "inst/extdata/dummy_species_rc.rds",
                       path_species = system.file("extdata", "dummy_species_rc.rds", package = "hibiiscus", mustWork = TRUE),
                       dataset_name = "custom - dummy_data")
  expect_true(inherits(tmp_db, 
                       "HibDatabase")
  )
})

