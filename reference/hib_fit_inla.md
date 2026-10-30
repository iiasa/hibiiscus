# hib_fit_inla: Fit an INLA mixed-effects model

This function fits an INLA model using the supplied database, formula
and other function arguments tuning the model fitting process

## Usage

``` r
hib_fit_inla(
  obj,
  response_fixed_form,
  random_covar,
  family = "gaussian",
  fit_settings,
  incl_spatial_term = FALSE,
  input_mesh = NULL,
  mesh_settings = list(mesh_dim = "2d", max.edge = c(75, 100), cutoff = 1e-12, offset =
    NULL),
  scale_specs = list(covars_to_scale = NULL, scale_func = NULL),
  relevel_specs = c(NULL),
  prior_list
)
```

## Arguments

- obj:

  A
  [`HibDatabase`](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
  object.

- response_fixed_form:

  A `formula` with the response and fixed effects

- random_covar:

  A `formula` with random, non-spatial, 'iid' effects specified

- family:

  A `character` indicating the model family, e.g., "gaussian"

- fit_settings:

  In development

- incl_spatial_term:

  A `logical`, whether or not to include a spatial random effect

- input_mesh:

  Optionally provide an INLA `mesh` object

- mesh_settings:

  Optionally specify INLA mash settings, e.g., max.edge

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

- prior_list:

  In development

## Value

HibModelINLA

## Details

This function takes the supplied database object (and various other
arguments), fits the specified model and returns a HibModelINLA object.
Still in development. Weighted data integration currently not
implemented.

## References

- Havard Rue, Sara Martino, and Nicholas Chopin (2009), Approximate
  Bayesian Inference for Latent Gaussian Models Using Integrated Nested
  Laplace Approximations (with discussion), Journal of the Royal
  Statistical Society B, 71, 319-39

- H. Rue, A. Riebler, S. H. Sorbye, J. B. Illian, D. P. Simpson,
  and F. K. Lindgren. Bayesian computing with INLA: A review. Annual
  Reviews of Statistics and Its Applications, 4(March):395-421, 2017.
  URL http://arxiv.org/abs/1604.00860

## See also

[HibModelINLA](https://iiasa.github.io/hibiiscus/reference/HibModelINLA-Class.md)
[`inla`](https://rdrr.io/pkg/INLA/man/inla.html)

## Examples

``` r
if (FALSE) { # \dontrun{
  mod <- hib_fit_inla(obj, y ~ -1 + Intercept + x, ~ f(a, model = "iid"))
} # } 
```
