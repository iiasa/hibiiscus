# Testthat checks for hib_filter.R

test_that("Check hib_filter correctly filters database",{
  tmp_db <- load_hib(#"inst/extdata/dummy_database.rds"
                     system.file("extdata", "dummy_database.rds", package = "hibiiscus", mustWork = TRUE))
  # -Select a random country, and exclude it
  tmp_db <- hib_filter_sites(tmp_db,
                             Country != tmp_db$sites$Country[1])
  # check the filter has worked
  expect_true(inherits(tmp_db, 
                       "Database")
  )
  expect_length(tmp_db$processing_flags$filtering_steps, 1)
  
  # -Select a random species, and exclude it
  tmp_db <- hib_filter_species(tmp_db,
                               Best_guess_binomial != tmp_db$species$Best_guess_binomial[1])
  # check the filter works
  expect_true(inherits(tmp_db, 
                       "Database")
  )
  expect_length(tmp_db$processing_flags$filtering_steps, 2)
  
})

