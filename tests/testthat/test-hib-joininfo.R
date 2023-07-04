# Testthat checks for hib_joininfo.R


tmp_db <- create_hib(path_species = "../raw_data/millard_et_al_2021/PREDICTS_pollinators_8_exp.rds",
                     version = "custom")
tmp_db <- hib_filter_sites(tmp_db, Country == "United Kingdom")

tmp <- "hello"

# %%%%%%%%%%
# hib_join_sites
# %%%%%%%%%%
# - 

tmp_db <- hib_join_sites(tmp_db,
                         data.frame(Predominant_land_use = c("Primary vegetation",                       
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
                         )

# - function works correctly
test_that("hib_join_sites returns a hib database", {
  expect_true(inherits(tmp_db, 
                       "Database")
  )
})


# - function fails if no shared columns
test_that("hib_join_sites returns an error if no shared columns", {
  expect_error(hib_join_sites(tmp_db,
                              data.frame(Not_a_column = c("Primary vegetation",                       
                                                                  "Mature secondary vegetation",
                                                                  "Young secondary vegetation",    
                                                                  "Secondary vegetation (indeterminate age)", 
                                                                  "Pasture",
                                                                  "Cropland", 
                                                                  "Urban"),
                                         Not_a_land_use = c("Natural",
                                                                 "Natural",
                                                                 "Semi-natural",
                                                                 "Semi-natural",
                                                                 "Human",
                                                                 "Human",
                                                                 "Human"))),
               regexp = "Columns in site_info are not compatible with the current site data"
  )
})

# - function fails if obj is not a database
test_that("hib_join_sites fails if obj is not a database", {
  expect_error(hib_join_sites(tmp,
                              data.frame(Predominant_land_use = c("Primary vegetation",                       
                                                                  "Mature secondary vegetation",
                                                                  "Young secondary vegetation",    
                                                                  "Secondary vegetation (indeterminate age)", 
                                                                  "Pasture",
                                                                  "Cropland", 
                                                                  "Urban"),
                                         New_column = c("Natural",
                                                                 "Natural",
                                                                 "Semi-natural",
                                                                 "Semi-natural",
                                                                 "Human",
                                                                 "Human",
                                                                 "Human"))),
               regexp = "obj does not inherit from class Database"
  )
})

# - function fails if site_info is not a dataframe
test_that("hib_join_sites fails if site_info is not a data frame", {
  expect_error(hib_join_sites(tmp_db,
                              "Not a dataframe"),
               regexp = "site_info is not a data frame"
  )
})


# %%%%%%%%%%
# hib_join_species
# %%%%%%%%%%

tmp_db <- create_hib(path_species = "../raw_data/millard_et_al_2021/PREDICTS_pollinators_8_exp.rds",
                     version = "custom")
tmp_db <- hib_filter_sites(tmp_db, Country == "United Kingdom")


# %%%
tmp_db <- hib_join_species(tmp_db,
                         data.frame(Predominant_land_use = c("Primary vegetation",                       
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
)

# - function works correctly
test_that("hib_join_species returns a hib database", {
  expect_true(inherits(tmp_db, 
                       "Database")
  )
})


# - function fails if no shared columns
test_that("hib_join_species returns an error if no shared columns", {
  expect_error(hib_join_species(tmp_db,
                              data.frame(Not_a_column = c("Primary vegetation",                       
                                                          "Mature secondary vegetation",
                                                          "Young secondary vegetation",    
                                                          "Secondary vegetation (indeterminate age)", 
                                                          "Pasture",
                                                          "Cropland", 
                                                          "Urban"),
                                         Not_a_land_use = c("Natural",
                                                            "Natural",
                                                            "Semi-natural",
                                                            "Semi-natural",
                                                            "Human",
                                                            "Human",
                                                            "Human"))),
               regexp = "Columns in spp_info are not compatible with the current species data"
  )
})

# - function fails if obj is not a database
test_that("hib_join_species fails if obj is not a database", {
  expect_error(hib_join_species(tmp,
                              data.frame(Predominant_land_use = c("Primary vegetation",                       
                                                                  "Mature secondary vegetation",
                                                                  "Young secondary vegetation",    
                                                                  "Secondary vegetation (indeterminate age)", 
                                                                  "Pasture",
                                                                  "Cropland", 
                                                                  "Urban"),
                                         New_column = c("Natural",
                                                        "Natural",
                                                        "Semi-natural",
                                                        "Semi-natural",
                                                        "Human",
                                                        "Human",
                                                        "Human"))),
               regexp = "obj does not inherit from class Database"
  )
})

# - function fails if site_info is not a dataframe
test_that("hib_join_species fails if site_info is not a data frame", {
  expect_error(hib_join_species(tmp_db,
                              "Not a dataframe"),
               regexp = "spp_info is not a data frame"
  )
})


