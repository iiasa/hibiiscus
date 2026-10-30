# Tests for hib_saveload.R


# %%%%%%%%%%
# hib_save
# %%%%%%%%%%

test_that("Check that hib_save correctly saves data.", {
  tmp_db <- hib_create_db(#path_species = "inst/extdata/dummy_species_rc.rds",
                       path_species = system.file("extdata", "dummy_species_rc.rds", package = "hibiiscus", mustWork = TRUE),
                       dataset_name = "custom - dummy_data")
  # - check it fails if not path provided
  expect_error(hib_save(tmp_db),
               regexp = "argument \"path\" is missing, with no default")
  # - fails for numeric
  expect_error(hib_save(tmp_db, 45),
               regexp = "path is not a character vector")
  
  # What about a non-hib object?!! - need to add
  
  # - And check we can actually save to file too
  # rm file if already exists??
  # save_hib(tmp_db,
  #          "inst/extdata/dummy_database.rds")
  # # check file exists now
  # expect_true(file.exists("inst/extdata/dummy_database.rds"))
})

test_that("hib_save and hib_load round-trip a single-raster AOA", {
  
  testthat::skip_if_not_installed("terra")
  
  projection <- terra::rast(nrows = 2, ncols = 2,
                            xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                            vals = 1:4)
  aoa <- terra::rast(nrows = 2, ncols = 2,
                     xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                     vals = c(0, 1, 1, 0))
  model <- structure(list(), class = "HibModelLME4")
  projection_obj <- HibProjection$new(
    obj = model,
    proj = projection,
    backtransform = NULL,
    rescale_specs = NULL,
    lme4_specs = NULL,
    inla_specs = NULL,
    brms_specs = NULL,
    aoa_res = list(proj_aoa = aoa)
  )

  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  save_path <- file.path(tmp_dir, "projection.rds")
  hib_save(projection_obj, save_path)
  loaded_obj <- hib_load(save_path)

  expect_s4_class(loaded_obj$projection, "SpatRaster")
  expect_s4_class(loaded_obj$aoa$proj_aoa, "SpatRaster")
  expect_equal(terra::values(loaded_obj$projection), terra::values(projection))
  expect_equal(terra::values(loaded_obj$aoa$proj_aoa), terra::values(aoa))

  unlink(tmp_dir, recursive = TRUE)
})


# %%%%%%%%%%
# hib_load
# %%%%%%%%%%

test_that("Check that hib_load correctly loads a pre-made hibiiscus database.", {
  # - fails if no path
  expect_error(hib_load(),
               regexp = "argument \"path\" is missing, with no default")
  # - fails for numeric
  expect_error(hib_load(45),
               regexp = "path is not a character vector")
  # - fails if not a hib database
  expect_error(hib_load(#"inst/extdata/dummy_species_rc.rds"
                        system.file("extdata", "dummy_species_rc.rds", package = "hibiiscus", mustWork = TRUE)),
               regexp = "is not TRUE")
  # - make sure it works for a hib db
  # tmp_db <- load_hib(#"inst/extdata/dummy_database.rds"
  #                    system.file("extdata", "dummy_database.rds", package = "hibiiscus", mustWork = TRUE))
  # expect_true(inherits(tmp_db, 
  #                      "HibDatabase"))
})  
