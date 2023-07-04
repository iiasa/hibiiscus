# Testthat checks for hib_lu_legend.R

# TODO: put millard data in inst/extdata/millard_et_al_2021_PREDICTS_pollinators_8_exp.rds

tmp_db <- create_hib(path_species = "../raw_data/millard_et_al_2021/PREDICTS_pollinators_8_exp.rds",
                     version = "custom")
tmp_db <- hib_filter_sites(tmp_db, Country == "United Kingdom")


alt_lu_df <- data.frame(Predominant_land_use = c("Primary vegetation",                       
                                                 "Mature secondary vegetation",
                                                 "Young secondary vegetation",    
                                                 "Secondary vegetation (indeterminate age)", 
                                                 "Pasture",
                                                 "Cropland", 
                                                 "Urban"),
                        Simplified_land_use = c("Natural",
                                                "Natural",
                                                "Semi-natural",
                                                "Semi-natural",
                                                "Human",
                                                "Human",
                                                "Human"))
alt_lu_df <- tmp_db$sites[,c("SSBS", "Predominant_land_use")] %>%
  dplyr::left_join(alt_lu_df) %>%
  dplyr::select(-Predominant_land_use) %>%
  dplyr::rename(Predominant_land_use = Simplified_land_use) %>%
  dplyr::mutate(Use_intensity = "",
                LUI = Predominant_land_use)


tmp_leg <- create_lu_legend(alt_lu_df)


# %%%%%%%%%%
# create_lu_legend
# %%%%%%%%%%
# - if input is a bad path, fail
test_that("create_lu_legend returns error if input is a bad file path", {
  expect_error(create_lu_legend("../not_a_path.rds"),
               regexp = "Path '../not_a_path.rds' does not exist"
  )
})

# - if input is not a df, fail
test_that("create_lu_legend returns error if input is not a df or file path", {
  expect_error(create_lu_legend(c(1,0,1)),
               regexp = "input is not a character vector or input is not a data frame"
  )
})

# - check colnames of output$lu_df - Predominant_land_use,  Use_intensity, LUI
test_that("create_lu_legend returns error if input doesn't contain correct LU columns", {
  expect_error(create_lu_legend(data.frame(SSBS = c("Example 1"),
                                           Predominant_land_use = c("Primary vegetation"))),
               regexp = "lu_df does not have all of these name"
  )
})

# - correct functioning
test_that("create_lu_legend returns an lu_legend", {
  expect_true(inherits(tmp_leg, 
                       "lu_legend")
  )
})


# %%%%%%%%%%
# save_legend
# %%%%%%%%%%

# - check it fails if not path provided
test_that("save_legend returns an error if no file path is provided", {
  expect_error(save_legend(tmp_leg),
               regexp = "argument \"output_path\" is missing, with no default")
})

# - fails for numeric
test_that("save_legend returns an error if file path is numeric", {
  expect_error(save_legend(tmp_leg, 45),
               regexp = "output_path is not a character vector")
})

# TODO: check that this also works...
# write to a temp file?? - write to "inst/extdata/example_lu_legend.rds"


# %%%%%%%%%%
# load_lu_legend
# %%%%%%%%%%

# - check it fails if path is wrong
test_that("load_lu_legend returns an error if file path is invlaid", {
  expect_error(load_lu_legend("../not_an_lu_legend.rds"),
               regexp = "Path '../not_an_lu_legend.rds' does not exist")
})

# - check it fails if file is not an lu legend
test_that("load_lu_legend returns an error if file is not a valid lu_legend", {
  expect_error(load_lu_legend("../processed_data/primary_vegetation_reclass.rds"),
               regexp = "leg does not inherit from class lu_legend")
})

# TODO: check that this also works...
# load from "inst/extdata/example_lu_legend.rds"

# %%%%%%%%%%
# set_lu_legend
# %%%%%%%%%%

# - only works if obj is a Database and leg is a lu_legend
test_that("set_lu_legend returns an error if obj is not a Database", {
  expect_error(set_lu_legend("hello",
                             tmp_leg),
               regexp = "obj does not inherit from class Database")
})

test_that("set_lu_legend returns an error if obj is not a Database", {
  expect_error(set_lu_legend(tmp_db,
                             c(1,0,1,0)),
               regexp = "leg does not inherit from class lu_legend")
})


tmp_db <- set_lu_legend(tmp_db,
                        tmp_leg)


# - correct functionality
test_that("set_lu_legend returns a hib database", {
  expect_true(inherits(tmp_db, 
                       "Database")
  )
})
