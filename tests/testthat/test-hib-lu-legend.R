# Testthat checks for hib_lu_legend.R

# %%%%%%%%%%
# create_lu_legend
# %%%%%%%%%%

test_that("Check that create_lu_legend works as intended.", {
  # - if input is a bad path, fail  
  expect_error(hib_create_lu_legend("data/not_a_path.rds"),
               regexp = "Path 'data/not_a_path.rds' does not exist"
  )
  # - if input is not a df, fail
  expect_error(hib_create_lu_legend(c(1,0,1)),
               regexp = "input is not a character vector or input is not a data frame"
  )
  # - check colnames of output$lu_df - Predominant_land_use,  Use_intensity, LUI
  expect_error(hib_create_lu_legend(data.frame(SSBS = c("Example 1"),
                                           Predominant_land_use = c("Primary vegetation"))),
               regexp = "lu_df does not have all of these name"
  )
  
  # - correct functioning
  alt_lu_df <- data.frame(SSBS = c("SSBS_1", "SSBS_2"),
                          Predominant_land_use = c("Natural/Semi-natural",                       
                                                   "Human"),
                          Use_intensity = c("Minimal_use",
                                            "Light_use")) |>
    dplyr::mutate(LUI = paste(Predominant_land_use, Use_intensity, 
                              sep = "_"))
  
  tmp_leg <- hib_create_lu_legend(alt_lu_df)
  expect_true(inherits(tmp_leg, 
                       "HibLULegend")
  )
  
})


# # %%%%%%%%%%
# # save_legend
# # %%%%%%%%%%

# test_that("Check that save_legend works correctly.", {
#   alt_lu_df <- data.frame(SSBS = c("SSBS_1", "SSBS_2"),
#                           Predominant_land_use = c("Natural/Semi-natural",                       
#                                                    "Human"),
#                           Use_intensity = c("Minimal_use",
#                                             "Light_use")) |>
#     dplyr::mutate(LUI = paste(Predominant_land_use, Use_intensity, 
#                               sep = "_"))
  
#   tmp_leg <- create_lu_legend(alt_lu_df)
  
#   # - check it fails if not path provided
#   expect_error(save_legend(tmp_leg),
#                regexp = "argument \"output_path\" is missing, with no default")
#   # - fails for numeric
#   expect_error(save_legend(tmp_leg, 45),
#                regexp = "output_path is not a character vector")
  
#   # TODO: check that this works...
#   # write to a temp file?? - write to "inst/extdata/example_lu_legend.rds"
  
# })


# # %%%%%%%%%%
# # load_lu_legend
# # %%%%%%%%%%

# test_that("Check load_lu_legend works as intended.", {
#   # - check it fails if path is wrong
#   expect_error(load_lu_legend("data/not_an_lu_legend.rds"),
#                regexp = "Path 'data/not_an_lu_legend.rds' does not exist")
#   # - check it fails if file is not an lu legend ...
#   expect_error(load_lu_legend(#"data/dummy_species_rc.rds"
#                               system.file("extdata", "dummy_species_rc.rds", package = "hibiiscus", mustWork = TRUE)),
#                regexp = "leg does not inherit from class HibLULegend")
#   # TODO: check that this also works...
#   # load from "inst/extdata/example_lu_legend.rds"
# })


# %%%%%%%%%%
# set_lu_legend
# %%%%%%%%%%
test_that("Check set_lu_legend works correctly.", {
  tmp_db <- hib_create_db(#path_species = "inst/extdata/dummy_species_rc.rds",
                       path_species = system.file("extdata", "dummy_species_rc.rds", package = "hibiiscus", mustWork = TRUE),
                       version = "custom",
                       custom_version_name = "dummy_data")
  # tmp_db <- load_hib(#"data/dummy_database.rds"
  #                    system.file("extdata", "dummy_database.rds", package = "hibiiscus", mustWork = TRUE))
  
  alt_lu_df <- tmp_db$sites[,c("SSBS", "Predominant_land_use", "Use_intensity")] |>
    dplyr::mutate(Predominant_land_use = dplyr::case_when(Predominant_land_use == "Primary" ~ "Natural/Semi-natural",
                                                          Predominant_land_use == "Secondary" ~ "Natural/Semi-natural",
                                                          Predominant_land_use == "Cropland" ~ "Human",
                                                          Predominant_land_use == "Urban" ~ "Human",
                                                          T ~ Predominant_land_use
                                                          ),
                  Use_intensity = dplyr::case_when(Use_intensity == "Minimal_use" ~ "Low",
                                                   Use_intensity == "Light_use" ~ "Low",
                                                   Use_intensity == "Intense_use" ~ "High",
                                                   T ~ Use_intensity 
                                                   ),
                  LUI = paste(Predominant_land_use, Use_intensity, 
                              sep = "_")) 
  
  tmp_leg <- hib_create_lu_legend(alt_lu_df)
  
  # - only works if obj is a Database and leg is a lu_legend
  expect_error(hib_set_lu_legend("hello",
                             tmp_leg),
               regexp = "obj does not inherit from class HibDatabase")
  expect_error(hib_set_lu_legend(tmp_db,
                             c(1,0,1,0)),
               regexp = "leg does not inherit from class HibLULegend")
  
  
  
  tmp_db <- hib_set_lu_legend(tmp_db,
                          tmp_leg)
  
  # - correct functionality
  expect_true(inherits(tmp_db, 
                       "HibDatabase"))
    
})


