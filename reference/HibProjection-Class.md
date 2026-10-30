# HibProjection

This is the template object for a an
[R6](https://r6.r-lib.org/reference/R6Class.html) class containing the
results of projecting a fitted HibModel based on data of interest. It
contains functions for summarising and visualising the results.

## Details

Template object for projection results

## Public fields

- `last_update`:

  (`character`)  
  Identifies the most recent update to the object (YYYY-MM-DD).

- `model`:

  (`R6`)  
  The model used to project, one of
  [HibModelLME4](https://iiasa.github.io/hibiiscus/reference/HibModelLME4-Class.md),
  [HibModelINLA](https://iiasa.github.io/hibiiscus/reference/HibModelINLA-Class.md),
  [HibModelBRMS](https://iiasa.github.io/hibiiscus/reference/HibModelBRMS-Class.md).

- `projection`:

  (`data.frame`, `SpatRaster`, or `list`)  
  The results of model projection.

- `projection_specs`:

  (`list`)  
  List of specifications used for modelling.

- `aoa`:

  (`list`)  
  Results of Area-of-Applicability calculations.

## Methods

### Public methods

- [`HibProjection$new()`](#method-HibProjection-initialize)

- [`HibProjection$print()`](#method-HibProjection-print)

- [`HibProjection$plot()`](#method-HibProjection-plot)

- [`HibProjection$clone()`](#method-HibProjection-clone)

------------------------------------------------------------------------

### `HibProjection$new()`

Creates a new instance of this
[R6](https://r6.r-lib.org/reference/R6Class.html) class. Note, this
object is constructed via the
[hib_project_model](https://iiasa.github.io/hibiiscus/reference/hib_project_model.md)
function found within this package.

#### Usage

    HibProjection$new(
      obj,
      proj,
      backtransform,
      rescale_specs,
      lme4_specs,
      inla_specs,
      brms_specs,
      aoa_res
    )

#### Arguments

- `obj`:

  (`HibModelXXXX`)  
  The model object used in cross-validation.

- `proj`:

  (`data.frame`, `SpatRaster`, or `list`)  
  The results of model projection.

- `backtransform`:

  (`function` or `NULL`)  
  Describes the backtransformation specified when making a model
  projection.

- `rescale_specs`:

  (`list`)  
  List of specifications used for rescaling predictions.

- `lme4_specs`:

  (`list`)  
  List of specifications used for lme4 models.

- `inla_specs`:

  (`list`)  
  List of specifications used for inla models.

- `brms_specs`:

  (`list`)  
  List of specifications used for brms models.

- `aoa_res`:

  (`list`, or `NULL`)  
  List of Area-of-Applicability results.

#### Returns

(`HibProjection`)  

------------------------------------------------------------------------

### `HibProjection$print()`

Base printout function presents core information about the object.

#### Usage

    HibProjection$print()

#### Returns

(`NULL`).

------------------------------------------------------------------------

### `HibProjection$plot()`

Function to generate figures summarising the projection results.

#### Usage

    HibProjection$plot(
      type = "response-map",
      response_var = NULL,
      uncertainty_var = NULL,
      aoa_mask = FALSE,
      aoa_var = NULL
    )

#### Arguments

- `type`:

  ([`character()`](https://rdrr.io/r/base/character.html)) Options are:
  "response-map", "uncertainty-map", "bivariate-map", "response-distr",
  "uncertainty-distr", "bivariate-distr".

- `response_var`:

  (`character`, or `NULL`)

- `uncertainty_var`:

  (`character`, or `NULL`)

- `aoa_mask`:

  (`logical`)

- `aoa_var`:

  (`character`, or `NULL`)

#### Returns

(`ggplot`).

------------------------------------------------------------------------

### `HibProjection$clone()`

The objects of this class are cloneable with this method.

#### Usage

    HibProjection$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
