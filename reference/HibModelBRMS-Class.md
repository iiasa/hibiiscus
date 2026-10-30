# HibModelBRMS

This is the template object for an
[R6](https://r6.r-lib.org/reference/R6Class.html) `brms` model class. It
contains functions for summarising and visualising the model results.

## Details

Template object for brms model

## Public fields

- `dataset_name`:

  (`character`)  
  Identifies the dataset used.

- `last_update`:

  (`character`)  
  Identifies the most recent update to the object (YYYY-MM-DD).

- `data`:

  (`data.frame`)  
  The data used in model fitting.

- `model`:

  (`brmsfit`)  
  The brms model object.

- `extra_args`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  Identifies addition arguments supplied during model fitting (found in
  '...' when calling `hib_fit_brms`).

- `fixedef_coefficients`:

  (`data.frame`)  
  The fixed effects coefficients.

- `metrics`:

  (`list`)  
  List containing model performance metrics.

- `processing_flags`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  Identifies the processing steps used for the HibDatabase object that
  was used for model fitting.

- `site_locations`:

  (`data.frame`)  
  Contains columns for "SSBS", "Longitude" ,"Latitude"

## Methods

### Public methods

- [`HibModelBRMS$new()`](#method-HibModelBRMS-initialize)

- [`HibModelBRMS$print()`](#method-HibModelBRMS-print)

- [`HibModelBRMS$plot()`](#method-HibModelBRMS-plot)

- [`HibModelBRMS$check_model()`](#method-HibModelBRMS-check_model)

- [`HibModelBRMS$clone()`](#method-HibModelBRMS-clone)

------------------------------------------------------------------------

### `HibModelBRMS$new()`

Creates a new instance of this
[R6](https://r6.r-lib.org/reference/R6Class.html) class. Note, this
object is constructed via the
[hib_fit_brms](https://iiasa.github.io/hibiiscus/reference/hib_fit_brms.md)
function found within this package.

#### Usage

    HibModelBRMS$new(
      db_obj,
      mod_df,
      brms_obj,
      fixef_df,
      metr_ls,
      n_spp,
      n_taxa,
      n_sites,
      extra_args
    )

#### Arguments

- `db_obj`:

  (`HibDatabase`, or `list`)  
  The HibDatabase object/s used to fit the model.

- `mod_df`:

  (`data.frame`)  
  The data.frame of the data used in model fitting.

- `brms_obj`:

  (`brmsfit`)  
  The brms model object.

- `fixef_df`:

  (`data.frame`)  
  The fixed effects from the fitted model.

- `metr_ls`:

  (`list`)  
  The list of model performance metrics.

- `n_spp`:

  (`integer` or `list`)  
  Number of species contributing to the fitted model.

- `n_taxa`:

  (`integer` or `list`)  
  Number of orignal taxa names contributing to the fitted model.

- `n_sites`:

  (`integer` or `list`)  
  Number of sites contributing to the fitted model.

- `extra_args`:

  (`list`)  
  Additional arguments specified during model fitting.

  @return (`HibModelBRMS`)

------------------------------------------------------------------------

### `HibModelBRMS$print()`

Base printout function presents core information about the object.

#### Usage

    HibModelBRMS$print()

#### Returns

(`NULL`).

------------------------------------------------------------------------

### `HibModelBRMS$plot()`

Function to generate figures summarising the model.

#### Usage

    HibModelBRMS$plot(type = "fixed-effects")

#### Arguments

- `type`:

  ([`character()`](https://rdrr.io/r/base/character.html)) Options are:
  "fixed-effects", "random-effects", "map-sites".

#### Returns

(`ggplot`).

------------------------------------------------------------------------

### `HibModelBRMS$check_model()`

Function to assess the model using the performance package.

#### Usage

    HibModelBRMS$check_model()

#### Returns

(`NULL`).

------------------------------------------------------------------------

### `HibModelBRMS$clone()`

The objects of this class are cloneable with this method.

#### Usage

    HibModelBRMS$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
