# Tests for createload_hib.R

# %%%%%%%%%%
# create_hib 
# %%%%%%%%%%


test_that("Check that create_hib correclty loads data.", {
  # - check that passing an incorrect file path returns an error
  expect_error(create_hib("inst/extdata/not_a_file.rds"),
               regexp = "does not exist")
  # - numeric path returns an error
  expect_error(create_hib(45),
               regexp = "not a character vector")
  # - ensure that incorrect version returns error
  expect_error(create_hib(#"inst/extdata/dummy_species_rc.rds",
                          path_species = system.file("extdata", "dummy_species_rc.rds", package = "hibiiscus", mustWork = TRUE),
                          version = "dummy_data"),
               regexp = "'arg' should be one of \"2016\", \"2022\", \"custom\"")
  
  # - check that only if correct files, and version are included will it work 
  tmp_db <- create_hib(#path_species = "inst/extdata/dummy_species_rc.rds",
                       path_species = system.file("extdata", "dummy_species_rc.rds", package = "hibiiscus", mustWork = TRUE),
                       version = "custom",
                       custom_version_name = "dummy_data")
  expect_true(inherits(tmp_db, 
                       "HibDatabase")
  )
})


# %%%%%%%%%%
# save_hib
# %%%%%%%%%%

test_that("Check that save_hib correctly saves data.", {
  tmp_db <- create_hib(#path_species = "inst/extdata/dummy_species_rc.rds",
                       path_species = system.file("extdata", "dummy_species_rc.rds", package = "hibiiscus", mustWork = TRUE),
                       version = "custom",
                       custom_version_name = "dummy_data")
  # - check it fails if not path provided
  expect_error(save_hib(tmp_db),
               regexp = "argument \"output_path\" is missing, with no default")
  # - fails for numeric
  expect_error(save_hib(tmp_db, 45),
               regexp = "output_path is not a character vector")
  
  # What about a non-hib object?!! - need to add
  
  # - And check we can actually save to file too
  # rm file if already exists??
  # save_hib(tmp_db,
  #          "inst/extdata/dummy_database.rds")
  # # check file exists now
  # expect_true(file.exists("inst/extdata/dummy_database.rds"))
})


# %%%%%%%%%%
# load_hib
# %%%%%%%%%%

test_that("Check that load_hib correctly loads a pre-made hibiiscus database.", {
  # - fails if no path
  expect_error(load_hib(),
               regexp = "argument \"path_input\" is missing, with no default")
  # - fails for numeric
  expect_error(load_hib(45),
               regexp = "path_input is not a character vector")
  # - fails if not a hib database
  expect_error(load_hib(#"inst/extdata/dummy_species_rc.rds"
                        system.file("extdata", "dummy_species_rc.rds", package = "hibiiscus", mustWork = TRUE)),
               regexp = "obj does not inherit from class HibDatabase")
  # - make sure it works for a hib db
  # tmp_db <- load_hib(#"inst/extdata/dummy_database.rds"
  #                    system.file("extdata", "dummy_database.rds", package = "hibiiscus", mustWork = TRUE))
  # expect_true(inherits(tmp_db, 
  #                      "HibDatabase"))
})  



