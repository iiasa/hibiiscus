# HibCV

This is the template object for a an
[R6](https://r6.r-lib.org/reference/R6Class.html) class containing the
results of cross-validation. It contains functions for summarising and
visualising the results.

## Details

Template object for cross-validation results

## Public fields

- `last_update`:

  (`character`)  
  Identifies the most recent update to the object (YYYY-MM-DD).

- `model_type`:

  (`character`)  
  Identifies the type of model; one of "HibModelLME4", "HibModelINLA",
  "HibModelBRMS".

- `original_model_obj`:

  (`list`)  
  Contains a list of core model information

- `cv_results`:

  (`list`)  
  Contains the list of results from cross-validation

- `insample_results`:

  (`list`)  
  Contains the list of results from in-sample validation

- `re_form`:

  (`list`)  
  A list of the re_form model structures used in evaluation.

- `use_original_mesh`:

  (`logical` or `NULL`)  
  Indicates use of original mesh for spatial INLA models.

## Methods

### Public methods

- [`HibCV$new()`](#method-HibCV-initialize)

- [`HibCV$print()`](#method-HibCV-print)

- [`HibCV$plot()`](#method-HibCV-plot)

- [`HibCV$clone()`](#method-HibCV-clone)

------------------------------------------------------------------------

### `HibCV$new()`

Creates a new instance of this
[R6](https://r6.r-lib.org/reference/R6Class.html) class. Note, this
object is constructed via the
[hib_cv_model](https://iiasa.github.io/hibiiscus/reference/hib_cv_model.md)
function found within this package.

#### Usage

    HibCV$new(mod_obj, cv_res, insamp_res, re_form, use_orig_mesh)

#### Arguments

- `mod_obj`:

  (`HibModelXXXX`)  
  The model object used in cross-validation.

- `cv_res`:

  (`list`)  
  The results of cross-validation.

- `insamp_res`:

  (`list`)  
  The results of ins-sample validation.

- `re_form`:

  (`list`)  
  The list of random effects structures used in model evaluation.

- `use_orig_mesh`:

  (`logical` or `NULL`)  
  Indicates use of original mesh for spatial INLA models

#### Returns

(`HibCV`)  

------------------------------------------------------------------------

### `HibCV$print()`

Base printout function presents core information about the object.

#### Usage

    HibCV$print()

#### Returns

(`NULL`).

------------------------------------------------------------------------

### `HibCV$plot()`

Function to generate figures summarising the cross- validation results.

#### Usage

    HibCV$plot(type = "metrics", metrics = c("rmse", "rsq"))

#### Arguments

- `type`:

  ([`character()`](https://rdrr.io/r/base/character.html)) Options are:
  "coefficients", "coefficients-summary", "obs-v-pred", "obs-v-error",
  "pred-v-error", "metrics", "map-sites".

- `metrics`:

  (`character`) Options are the lower-case abbreviations of calculated
  error metrics. - Only needed if type == "metrics".

#### Returns

(`ggplot`).

------------------------------------------------------------------------

### `HibCV$clone()`

The objects of this class are cloneable with this method.

#### Usage

    HibCV$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
