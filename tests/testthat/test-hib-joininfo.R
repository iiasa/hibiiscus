# Testthat checks for hib_joininfo.R


# %%%%%%%%%%
# hib_join_sites
# %%%%%%%%%%

test_that("Check that hib_join_sites works correctly.", {
  tmp_db <- load_hib(#"inst/extdata/dummy_database.rds"
                     system.file("extdata", "dummy_database.rds", package = "hibiiscus", mustWork = TRUE))
  # tmp_db$sites$Predominant_land_use |> unique()
  # Default levels of dummy data are: "Primary" "Secondary" "Cropland" "Urban"
  # Simplified levels: "Natural/Semi-natural" "Human"
  
  tmp_df <- data.frame(Predominant_land_use = c("Primary", "Secondary", "Cropland", "Urban"),
                   Simplified_land_use = c("Natural/Semi-natural",
                                           "Natural/Semi-natural",
                                           "Human",
                                           "Human"),
                   Another_new_column = c("Blahhh"))
  
  tmp_db <- hib_join_sites(tmp_db,
                           tmp_df)
  # tmp_db$processing_flags
  # -Successful join
  expect_true(inherits(tmp_db, 
                       "Database")
  )
  # -No compatible columns
  expect_error(hib_join_sites(tmp_db,
                              tmp_df |>
                                dplyr::select(-Another_new_column) |> 
                                dplyr::rename(Not_a_column = Predominant_land_use,
                                              Not_a_landuse = Simplified_land_use)),
               regexp = "Columns in site_info are not compatible with the current site data"
               )
  # -Not a database for obj
  expect_error(hib_join_sites("Hello",
                              tmp_df),
               regexp = "obj does not inherit from class Database"
  )
  # -Not a dataframe for site_info
  expect_error(hib_join_sites(tmp_db,
                              "Not a dataframe"),
               regexp = "site_info is not a data frame"
  )
})


# %%%%%%%%%%
# hib_join_species
# %%%%%%%%%%

test_that("Check that hib_join_species works correctly.", {
  tmp_db <- load_hib(#"inst/extdata/dummy_database.rds"
                     system.file("extdata", "dummy_database.rds", package = "hibiiscus", mustWork = TRUE))
  # tmp_db$sites$Predominant_land_use |> unique()
  # Default levels of dummy data are: "Primary" "Secondary" "Cropland" "Urban"
  # Simplified levels: "Natural/Semi-natural" "Human"
  
  tmp_df <- data.frame(Predominant_land_use = c("Primary", "Secondary", "Cropland", "Urban"),
                   Simplified_land_use = c("Natural/Semi-natural",
                                           "Natural/Semi-natural",
                                           "Human",
                                           "Human"))
  
  tmp_db <- hib_join_species(tmp_db,
                             tmp_df)
  # tmp_db$processing_flags
  # -Successful join
  expect_true(inherits(tmp_db, 
                       "Database")
  )
  # -No compatible columns
  expect_error(hib_join_species(tmp_db,
                                tmp_df |>
                                dplyr::rename(Not_a_column = Predominant_land_use,
                                              Not_a_landuse = Simplified_land_use)),
               regexp = "Columns in spp_info are not compatible with the current species data"
  )
  # -Not a database for obj
  expect_error(hib_join_species("Hello",
                                tmp_df),
               regexp = "obj does not inherit from class Database"
  )
  # -Not a dataframe for site_info
  expect_error(hib_join_species(tmp_db,
                              "Not a dataframe"),
               regexp = "spp_info is not a data frame"
  )
})


