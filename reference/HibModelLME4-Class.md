# HibModelLME4

This is the template object for an
[R6](https://r6.r-lib.org/reference/R6Class.html) `lme4` model class. It
contains functions for summarising and visualising the model results.

## Details

Template object for lme4 model

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

  (`merMod`)  
  The lme4 model object.

- `fixedef_coefficients`:

  (`data.frame`)  
  The fixed effects coefficients.

- `ci_type`:

  (`character`, or `NA`)  
  Identifies the types of confidence intervals calculated.

- `metrics`:

  (`numeric`)  
  Named vector containing model performance metrics.

- `processing_flags`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  Identifies the processing steps used for the HibDatabase object that
  was used for model fitting.

- `site_locations`:

  (`data.frame`)  
  Contains columns for "SSBS", "Longitude" ,"Latitude".

## Methods

### Public methods

- [`HibModelLME4$new()`](#method-HibModelLME4-initialize)

- [`HibModelLME4$print()`](#method-HibModelLME4-print)

- [`HibModelLME4$plot()`](#method-HibModelLME4-plot)

- [`HibModelLME4$check_model()`](#method-HibModelLME4-check_model)

- [`HibModelLME4$clone()`](#method-HibModelLME4-clone)

------------------------------------------------------------------------

### `HibModelLME4$new()`

Creates a new instance of this
[R6](https://r6.r-lib.org/reference/R6Class.html) class. Note, this
object is constructed via the
[hib_fit_lme4](https://iiasa.github.io/hibiiscus/reference/hib_fit_lme4.md)
function found within this package.

#### Usage

    HibModelLME4$new(
      db_obj,
      mod_df,
      lme4_obj,
      fixef_df,
      ci_type,
      metr_vec,
      n_spp,
      n_taxa,
      n_sites
    )

#### Arguments

- `db_obj`:

  (`HibDatabase`, or `list`)  
  The HibDatabase object/s used to fit the model.

- `mod_df`:

  (`data.frame`)  
  The data.frame of the data used in model fitting.

- `lme4_obj`:

  (`merMod`)  
  The lme4 model object.

- `fixef_df`:

  (`data.frame`)  
  The fixed effects from the fitted model.

- `ci_type`:

  (`character`, or `NA`)  
  The type of confidence intervals calculated.

- `metr_vec`:

  (`numeric`)  
  The named vector of model performance metrics.

- `n_spp`:

  (`integer` or `list`)  
  Number of species contributing to the fitted model.

- `n_taxa`:

  (`integer` or `list`)  
  Number of orignal taxa names contributing to the fitted model.

- `n_sites`:

  (`integer` or `list`)  
  Number of sites contributing to the fitted model.

  @return (`HibModelLME4`)

------------------------------------------------------------------------

### `HibModelLME4$print()`

Base printout function presents core information about the object.

#### Usage

    HibModelLME4$print()

#### Returns

(`NULL`).

------------------------------------------------------------------------

### `HibModelLME4$plot()`

Function to generate figures summarising the model.

#### Usage

    HibModelLME4$plot(type = "fixed-effects")

#### Arguments

- `type`:

  ([`character()`](https://rdrr.io/r/base/character.html)) Options are:
  "fixed-effects", "random-effects", "map-sites".

#### Returns

(`ggplot`).

------------------------------------------------------------------------

### `HibModelLME4$check_model()`

Function to assess the model using the performance package.

#### Usage

    HibModelLME4$check_model()

#### Returns

(`NULL`).

------------------------------------------------------------------------

### `HibModelLME4$clone()`

The objects of this class are cloneable with this method.

#### Usage

    HibModelLME4$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
