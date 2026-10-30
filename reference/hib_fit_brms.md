# hib_fit_brms: Fit a brms style mixed-effects model

This function fits a brms model using the supplied database, formula and
other function arguments tuning the model fitting process.

## Usage

``` r
hib_fit_brms(
  obj,
  form,
  family = stats::gaussian,
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
  object or a list of them.

- form:

  A `formula` in the style of `brms` (as in `lme4`) - including fixed
  and random effects.

- family:

  A `family` function, e.g. gaussian, poisson - as allowed by brms

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
  Note: if using this argument, the supplied formula should have a
  response of the form `y | weights(data_weight) ~`. The data_weight
  variable will be generated within the function call.

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
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html),
  e.g., prior, data2

## Value

`HibModelBRMS`

## Details

This function takes the supplied database object (and various other
arguments), fits the specified model and returns a HibModelBRMS object.

## References

- Paul-Christian Bürkner (2017). brms: An R Package for Bayesian
  Multilevel Models Using Stan. Journal of Statistical Software, 80(1),
  1-28. doi:10.18637/jss.v080.i01

- Paul-Christian Bürkner (2018). Advanced Bayesian Multilevel Modeling
  with the R Package brms. The R Journal, 10(1), 395-411.
  doi:10.32614/RJ-2018-017

- Paul-Christian Bürkner (2021). Bayesian Item Response Modeling in R
  with brms and Stan. Journal of Statistical Software, 100(5), 1-54.
  doi:10.18637/jss.v100.i05

## See also

[HibModelBRMS](https://iiasa.github.io/hibiiscus/reference/HibModelBRMS-Class.md)
[`brm`](https://paulbuerkner.com/brms/reference/brm.html)

## Examples

``` r
if (FALSE) { # \dontrun{
  mod <- hib_fit_brms(obj, y ~ x + ( 1 | a ), ...)
} # } 
```
