# Testthat checks for hib_filter.R


tmp_db <- create_hib(path_species = "../raw_data/millard_et_al_2021/PREDICTS_pollinators_8_exp.rds",
                     version = "custom")
tmp_db_filtered1 <- hib_filter_sites(tmp_db, Country == "United Kingdom")
tmp_db_filtered2 <- hib_filter_species(tmp_db, Country == "United Kingdom")

tmp <- "hello"
# %%%%%%%%%%
# hib_filter_sites
# %%%%%%%%%%
# - 
test_that("hib_filter_sites returns a hib database", {
  expect_true(inherits(tmp_db_filtered1, 
                       "Database")
  )
})


test_that("hib_filter_sites returns an error if input is not a Database", {
  expect_error(hib_filter_sites(tmp, Country == "United Kingdom"),
               regexp = "obj does not inherit from class Database"
  )
})




# %%%%%%%%%%
# hib_filter_species 
# %%%%%%%%%%
# - 
test_that("hib_filter_species returns a hib database", {
  expect_true(inherits(tmp_db_filtered2, 
                       "Database")
              )
})


test_that("hib_filter_species returns an error if input is not a Database", {
  expect_error(hib_filter_species(tmp, Country == "United Kingdom"),
               regexp = "obj does not inherit from class Database"
  )
})

