# HibModelINLA

This is the template object for an
[R6](https://r6.r-lib.org/reference/R6Class.html) `INLA` model class. It
contains functions for summarising and visualising the model results.

## Details

Template object for INLA model

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

  (`inla`)  
  The [`INLA::inla`](https://rdrr.io/pkg/INLA/man/inla.html) model
  object.

- `processing_flags`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  Identifies the processing steps used for the HibDatabase object that
  was used for model fitting.

- `fixedef_coefficients`:

  (`data.frame`)  
  The fixed effects coefficients.

- `user_formula`:

  (`formula` or `list`)  
  The model formula supplied by the user (includes the fixed effects and
  non-spatial random effects).

- `inla_fitted_formula`:

  (`formula`)  
  The fitted model formula (includes all effects).

- `inla_incl_spatial_term`:

  (`logical`)  
  TRUE or FALSE, does the model include a spatial term.

- `inla_stack`:

  (`inla.data.stack`)  
  The data.stack used in model fitting, an INLA object.

- `inla_mesh`:

  (`inla.mesh` or `NULL`)  
  If used, the mesh used in model fitting.

- `inla_mesh_settings`:

  (`list`)  
  A list of specifications for generating the mesh, if used in
  modelling.

- `inla_mesh_source`:

  (`character` or `NULL`)  
  A description of the mesh source - if not null, can be either
  "user_provided" or "mesh_settings".

- `site_locations`:

  (`data.frame`)  
  Contains columns for "SSBS", "Longitude" ,"Latitude".

## Methods

### Public methods

- [`HibModelINLA$new()`](#method-HibModelINLA-initialize)

- [`HibModelINLA$print()`](#method-HibModelINLA-print)

- [`HibModelINLA$plot()`](#method-HibModelINLA-plot)

- [`HibModelINLA$clone()`](#method-HibModelINLA-clone)

------------------------------------------------------------------------

### `HibModelINLA$new()`

Creates a new instance of this
[R6](https://r6.r-lib.org/reference/R6Class.html) class. Note, this
object is constructed via the
[hib_fit_inla](https://iiasa.github.io/hibiiscus/reference/hib_fit_inla.md)
function found within this package.

#### Usage

    HibModelINLA$new(
      db_obj,
      mod_df,
      inla_obj,
      fixef_df,
      user_formula,
      inla_fitted_formula,
      inla_incl_spatial_term,
      inla_stack,
      inla_mesh,
      inla_mesh_settings,
      inla_mesh_source,
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

- `inla_obj`:

  (`inla`)  
  The inla model object.

- `fixef_df`:

  (`data.frame`)  
  The fixed effects from the fitted model.

- `user_formula`:

  (`formula` or `list`)  
  The model formula supplied by the user.

- `inla_fitted_formula`:

  (`formula`)  
  The fitted model formula.

- `inla_incl_spatial_term`:

  (`logical`)  
  TRUE or FALSE, does the model include a spatial term.

- `inla_stack`:

  (`inla.data.stack`)  
  The data.stack used in model fitting, an INLA object.

- `inla_mesh`:

  (`inla.mesh` or `NULL`)  
  If used, the mesh used in model fitting.

- `inla_mesh_settings`:

  (`list`)  
  A list of specifications for generating the mesh.

- `inla_mesh_source`:

  (`character` or `NULL`)  
  A description of the mesh source.

- `n_spp`:

  (`integer`)  
  Number of species contributing to the fitted model.

- `n_taxa`:

  (`integer`)  
  Number of orignal taxa names contributing to the fitted model.

- `n_sites`:

  (`integer`)  
  Number of sites contributing to the fitted model.

  @return (`HibModelINLA`)

- `metr_ls`:

  (`list`)  
  The list of model performance metrics.

------------------------------------------------------------------------

### `HibModelINLA$print()`

Base printout function presents core information about the object.

#### Usage

    HibModelINLA$print()

#### Returns

(`NULL`).

------------------------------------------------------------------------

### `HibModelINLA$plot()`

Function to generate figures summarising the model.

#### Usage

    HibModelINLA$plot(type = "fixed-effects")

#### Arguments

- `type`:

  ([`character()`](https://rdrr.io/r/base/character.html)) Options are:
  "fixed-effects", "random-effects", "spatial-random-field",
  "map-sites".

#### Returns

(`ggplot`).

------------------------------------------------------------------------

### `HibModelINLA$clone()`

The objects of this class are cloneable with this method.

#### Usage

    HibModelINLA$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
