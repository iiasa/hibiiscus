# hib_cv_model: Cross-validate a fitted hibiiscus mixed-effects model

This function re-fits models on different subsets of the data, and
estimates a range of error metrics allowing for evaluation of model
predictive performance.

## Usage

``` r
hib_cv_model(
  mod_obj,
  db_obj,
  fold_type,
  custom_col = NULL,
  random_seed = NULL,
  re_form = list(NA),
  backtransform = NULL,
  use_original_mesh = NULL
)
```

## Arguments

- mod_obj:

  One of;
  [`HibModelLME4`](https://iiasa.github.io/hibiiscus/reference/HibModelLME4-Class.md),
  [`HibModelINLA`](https://iiasa.github.io/hibiiscus/reference/HibModelINLA-Class.md)
  or
  [`HibModelBRMS`](https://iiasa.github.io/hibiiscus/reference/HibModelBRMS-Class.md).

- db_obj:

  A
  [`HibDatabase`](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
  object, or, if the above model used multiple data sources, a `list` of
  `HibDatabase` objects.

- fold_type:

  A `character` indicating the type of cross-validation; one of:
  "Random-2", "Random-3", "Random-5", "Random-10", "Source_ID", "SS",
  "Country", "UN_subregion", "UN_region", "Biome", "Realm", "Custom".

- custom_col:

  A `character` indicating the column in the model data (`sites`
  data.frame) to use to define cross-validation folds, when
  `fold_type == "Custom"`.

- random_seed:

  A `numeric` indicating the random seed used to define random
  cross-validation folds - only used if 'Random-X' in `fold_type`.

- re_form:

  A `list` indicating the random effects variables to be included in the
  cross-validation/insample predicitons - if multiple items supplied,
  each option will be used in a separate cross-validated prediction.
  Values can be either `NA` (no random effects), `NULL` (all random
  effects), or a rhs `formula` specifying particular random effects to
  include.  
  Note: only test set data points with random-effect values present in
  the training data will have predictions if `re_form` is not `NA`.

- backtransform:

  A function to convert the model prediction back to the original scale
  of interest. The inverse of any transformation used in model fitting.

- use_original_mesh:

  A `logical` value indicating whether cross-validation of spatial
  `INLA` models should use the original mesh for all CV folds.
  Computationally this is slightly more efficient, but will impact CV
  performance metrics as 'withheld' sites will still exist in the mesh.

## Value

HibCV

## Details

This function allows for model evaluation using multiple
cross-validation strategies, and multiple random effects structures. It
returns a HibCV object, which includes the original model, insample and
CV performance for multiple evaluation metrics, model coefficients
estimated per cross-validation fold.

## See also

[HibCV](https://iiasa.github.io/hibiiscus/reference/HibCV-Class.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  cv_out <- hib_cv_model(mod_obj, db_obj, fold_type = "Random-10")
} # } 
```
