# hib_project_model: Project a hibiiscus mixed-effects model

This function projects a fitted HibModelXXX using the supplied model
object and data to be used in projection.

## Usage

``` r
hib_project_model(
  mod_obj,
  data,
  backtransform = NULL,
  rescale_specs = NULL,
  lme4_specs = NULL,
  inla_specs = NULL,
  brms_specs = NULL,
  calc_aoa = FALSE,
  cv_obj = NULL
)
```

## Arguments

- mod_obj:

  One of;
  [`HibModelLME4`](https://iiasa.github.io/hibiiscus/reference/HibModelLME4-Class.md),
  [`HibModelINLA`](https://iiasa.github.io/hibiiscus/reference/HibModelINLA-Class.md)
  or
  [`HibModelBRMS`](https://iiasa.github.io/hibiiscus/reference/HibModelBRMS-Class.md).

- data:

  A `data.frame`-like object, a `SpatRaster`, or a list of `SpatRaster`.
  Should contain all fixed-effect variables used in the fitted model,
  and any required random variables. If `data` is a `data.frame`, `x`
  and `y` columns should also be included.

- backtransform:

  A `function` that will backtransform the model response to the scale
  of the original variable, e.g.,
  [`exp()`](https://rdrr.io/r/base/Log.html) as the inverse of
  [`log()`](https://rdrr.io/r/base/Log.html). Default is `NULL`, i.e.,
  no backtransformation.

- rescale_specs:

  A named `list` (or `NULL`) containing details relevant to rescaling
  projections. Specifically, named items should consist of:

  - `rescale` (`logical`),

  - `baseline_covar` (`data.frame`, or NULL),

  - `baseline_re_form` (`formula`, `NA`, or `NULL`).  

  In combination, these arguments determine if predictions will be
  rescaled or not (`rescale`), the covariate values used to generate the
  baseline against which other predictions are rescaled
  (`baseline_covar`), and whether/how random effects should be
  incorporated in the generation of this baseline (`baseline_re_form`,
  see `re_form` in `lme4/inla/brms_specs` argument).  
  If `rescale == TRUE`, but no `baseline_covar` is provided, the
  intercept of the model will be used as the baseline value.

- lme4_specs:

  A named `list` (or `NULL`) containing details relevant to projecting a
  `HibModelLME4`.  
  Named items should consist of:

  - `bootstrap` (`logical`), indicating whether bootstrapped prediction
    intervals should be estimated. Can be slow for large models.

  - `n_boot` (`numeric`), how many bootstrap samples should be used to
    estimate prediction uncertainty.

  - `calc_se` (`logical`), specify whether or not a standard error based
    uncertainty estimate (via
    [`marginaleffects::predictions`](https://rdrr.io/pkg/marginaleffects/man/predictions.html))
    is calculated and used to generate prediction intervals.

  - `conf_level` (`numeric`) is used to specify the confidence levels of
    intervals generated via bootstrapping or the SE-based approach.
    Values should fall between 0 and 1, typically in the range 0.5 to
    0.95. Multiple values can be provided.

  - `re_form` (`formula`, `NULL`, or `NA`), as `re.form` in
    [`predict.merMod()`](https://rdrr.io/pkg/lme4/man/predict.merMod.html),
    used to specify random effects to condition on when predicting.

- inla_specs:

  A named `list` (or `NULL`) containing details relevant to projecting a
  `HibModelINLA`.  

  - `re_form` (`formula`, `NULL`, or `NA`), as `re.form` in
    [`predict.merMod()`](https://rdrr.io/pkg/lme4/man/predict.merMod.html),
    used to specify random effects to condition on when predicting.

  - `link_inverse` a `function` that is the inverse of the link used to
    fit the model. Used to backtranform the global intercept for use in
    rescaling predictions.

- brms_specs:

  A named `list` (or `NULL`) containing details relevant to projecting a
  `HibModelBRMS`. Named items should consist of:

  - `incl_obs_uncert` (`logical`), to determine whether or not
    observation-level uncertainty is included in predictions. Setting
    this parameter to `TRUE` will produce predictions using
    [`predict.brmsfit()`](https://paulbuerkner.com/brms/reference/predict.brmsfit.html)
    (larger uncertainty), setting to `FALSE` calls
    [`fitted.brmsfit()`](https://paulbuerkner.com/brms/reference/fitted.brmsfit.html)
    internally.

  - `conf_level` (`numeric`) is used to specify the confidence levels of
    intervals generated. Values should fall between 0 and 1, typically
    in the range 0.5 to 0.95. Internally converted to quantiles and
    passed to the `probs` argument of `predict.brmsfit()`/
    `fitted.brmsfit()`.

  - `re_form` (`formula`, `NULL`, or `NA`), as `re_formula` in
    `predict.brmsfit()`, used to specify random effects to condition on
    when predicting.

- calc_aoa:

  A `logical` to indicate whether or not the 'Area-of-Applicability'
  (AOA) should be estimated - following Meyer and Pebesma (2021).

- cv_obj:

  A
  [`HibCV`](https://iiasa.github.io/hibiiscus/reference/HibCV-Class.md)
  object, used in AOA estimation.

## Value

HibProjection

## Details

Given a HibModelXXX, new data and specification, this function
projects/predicts modelled responses for each 'unit' of the new data.
Calculation of the 'Area-of-Applicability' of the model within the new
data is an optional addition. Specification arguments allow for
customisation of the prediction type. If specifications are not
provided, defaults will be used. Projections will be returned in the
same format as the provided data. If data is provided as `SpatRaster`, a
land-use variable should be included as a single categorical layer, or
as multiple, proportional layers, named to align with the levels of the
categorical variable. Currently, only a single 'categorical' variable
can be provided in 'proportional cover' format in a `SpatRaster`.

## References

- Meyer, H., & Pebesma, E. (2021). Predicting into unknown space?
  Estimating the area of applicability of spatial prediction models.
  Methods in Ecology and Evolution, 12, 1620–1633.
  https://doi.org/10.1111/2041-210X.13650

- Meyer H, Milà C, Ludwig M, Linnenbrink J (2023). *CAST: 'caret'
  Applications for Spatial-Temporal Models*. R package,
  <https://CRAN.R-project.org/package=CAST>.

## See also

[`HibProjection`](https://iiasa.github.io/hibiiscus/reference/HibProjection-Class.md)
[`hib_calc_aoa`](https://iiasa.github.io/hibiiscus/reference/hib_calc_aoa.md)
[`predict.merMod`](https://rdrr.io/pkg/lme4/man/predict.merMod.html)
[`fitted.brmsfit`](https://paulbuerkner.com/brms/reference/fitted.brmsfit.html)
[`predict.brmsfit`](https://paulbuerkner.com/brms/reference/predict.brmsfit.html)

## Examples

``` r
if (FALSE) { # \dontrun{
  proj <- hib_mod_project(hib_model, proj_data)
} # } 
```
