# Tests for createload_hib.R

# %%%%%%%%%%
# create_hib 
# %%%%%%%%%%
# - check that passing an incorrect file path returns an error
test_that("create_hib returns an error if file does not exist", {
  expect_error(create_hib("../not_a_file.rds"),
               regexp = "does not exist")
})

# - numeric path returns an error
test_that("create_hib returns an error if file is not character", {
  expect_error(create_hib(45),
               regexp = "not a character vector")
})

# - ensure that incorrect version returns error
test_that("create_hib returns an error if version is not valid", {
  expect_error(create_hib("../raw_data/millard_et_al_2021/PREDICTS_pollinators_8_exp.rds",
                          version = "Millard_et_al_2021"),
               regexp = "'arg' should be one of \"2016\", \"2022\", \"custom\"")
})

# - return an error if file exists, version is correct, but file/s is not a PREDICTS-like database/not rds
test_that("create_hib returns an error if file is not RDS", {
  expect_error(create_hib("../raw_data/wwf_ecoregions/layerfiles_readme.txt",
                          version = "custom"),
               regexp = "does not have extension rds")
})

test_that("create_hib returns an error if file is not PREDICTS-like", {
  expect_error(create_hib("../processed_data/primary_vegetation_reclass.rds",
                          version = "custom"),
               regexp = "species does not have all of these name")
})

# - check that only if correct files, and version  are included will it work - conditional on class_database.
tmp_db <- create_hib(path_species = "../raw_data/millard_et_al_2021/PREDICTS_pollinators_8_exp.rds",
                     version = "custom")
test_that("create_hib returns a 'Database'", {
  expect_true(inherits(tmp_db, 
                       "Database")
              )
})


# %%%%%%%%%%
# save_hib
# %%%%%%%%%%
# - check it fails if not path provided
test_that("save_hib returns an error if no file path is provided", {
  expect_error(save_hib(tmp_db),
               regexp = "argument \"output_path\" is missing, with no default")
})

# - fails for numeric
test_that("save_hib returns an error if file path is numeric", {
  expect_error(save_hib(tmp_db, 45),
               regexp = "output_path is not a character vector")
})

# - fails if directory doesn't exist... - this needs to be fixed...

# %%%%%%%%%%
# load_hib
# %%%%%%%%%%
# - fails if no path
test_that("load_hib returns an error if no file path is provided", {
  expect_error(load_hib(),
               regexp = "argument \"path_input\" is missing, with no default")
})

# - fails for numeric
test_that("load_hib returns an error if file path is numeric", {
  expect_error(load_hib(45),
               regexp = "path_input is not a character vector")
})

# - fails if not a hib database
test_that("load_hib returns an error if file is not PREDICTS-like", {
  expect_error(load_hib("../processed_data/primary_vegetation_reclass.rds"),
               regexp = "obj does not inherit from class Database")
})

# - make sure it works for a hib db
tmp_hib <- load_hib("../processed_data/hib_millard_diversity_20230404.rds")
test_that("load_hib returns an error if file path is numeric", {
  expect_true(inherits(tmp_hib, 
                       "Database"))
})


