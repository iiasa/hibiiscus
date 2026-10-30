# hib_fit_lme4: Fit an lme4 style mixed-effects model

This function fits an lme4 model using the supplied database, formula
and other function arguments tuning the model fitting process.

## Usage

``` r
hib_fit_lme4(
  obj,
  form,
  family = stats::gaussian,
  REML = TRUE,
  calc_CIs = NULL,
  optimise_ranefs = FALSE,
  model_dredge = FALSE,
  model_average = FALSE,
  weight_specs = list(weights = NULL, weight_type = NULL),
  scale_specs = list(covars_to_scale = NULL, scale_func = NULL),
  relevel_specs = c(NULL),
  ...
)
```

## Arguments

- obj:

  A
  [`HibDatabase`](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
  object, or a list of databases.

- form:

  Model formula in the style of lme4 - including fixed and random
  effects.

- family:

  A `family` function, e.g. gaussian, poisson - as allowed by lme4.
  Dafault is 'gaussian'.

- REML:

  A `logical` - if using a gaussian family, should the model be fitted
  using restricted maximum likelihood (REML, T) or maximum
  likelihood (F) - as in
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html).

- calc_CIs:

  An indication of whether fixed effect confidence intervals (95%)
  should be calculated or not, and if so how - one of NULL, "Wald",
  "profile", "boot".

- optimise_ranefs:

  Not currently applicable - future dev

- model_dredge:

  Not currently applicable - future dev

- model_average:

  Not currently applicable - future dev

- weight_specs:

  A named `list` containing "weights" and "weight_type".  

  - `weights`: a `numeric` vector containing the relative weights to
    assign to each dataset/site, should sum to 1.  

  - `weight_type`: a `character`, specifying either "dataset" or
    "site" - i.e. should the weights argument assign relative weights to
    whole datasets, or for sites in each dataset.  

  Only applicable if a list of `HibDatabase`s are provided.  
  The length of 'weights' should be equal to the number of
  `HibDatabase`s.

- scale_specs:

  A named `list` containing "covars_to_scale" and "scale_func".

  - `covars_to_scale` is a `character` vector of numeric covariates that
    should be scaled prior to modelling.

  - `scale_func` is a function that will be used to scale the specified
    covariates - if `NULL`,
    [`base::scale()`](https://rdrr.io/r/base/scale.html) will be used.

- relevel_specs:

  A `vector` where each item is:
  `"<covar_to_relevel_name>" = "<base_level>"`.

- ...:

  Any other named arguments that can be passed to
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)/[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html),
  e.g., na.action

## Value

`HibModelLME4`

## Details

This function takes the supplied database object (and various other
arguments), fits the specified model and returns a HibModelLME4 object.

## References

- Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical
  Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.

## See also

[HibModelLME4](https://iiasa.github.io/hibiiscus/reference/HibModelLME4-Class.md)
[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html)
[`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html)

## Examples

``` r
if (FALSE) { # \dontrun{
  mod <- hib_fit_lme4(obj, y ~ x + ( 1 | a ), ...)
} # } 
```
