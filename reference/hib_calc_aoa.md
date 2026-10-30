# hib_calc_aoa: Calculate the Area-of-Applicability.

Calculate the Area-of-Applicability (AOA) given a model,
cross-validation output, and the data being used to make
projection/prediction.

## Usage

``` r
hib_calc_aoa(mod_obj, cv_obj, proj_data)
```

## Arguments

- mod_obj:

  One of;
  [`HibModelLME4`](https://iiasa.github.io/hibiiscus/reference/HibModelLME4-Class.md),
  [`HibModelINLA`](https://iiasa.github.io/hibiiscus/reference/HibModelINLA-Class.md)
  or
  [`HibModelBRMS`](https://iiasa.github.io/hibiiscus/reference/HibModelBRMS-Class.md).

- cv_obj:

  A
  [`HibCV`](https://iiasa.github.io/hibiiscus/reference/HibCV-Class.md)
  object, based on the above `mod_obj`.

- proj_data:

  A `data.frame`, `SpatRaster` or list of `SpatRaster`, containing
  covariates for projecting/predicting the model. Predicted responses
  should not be included.

## Value

list

## Details

This function estimates an AOA, for each of the cross-validation
approaches supplied. We follow the methods of Meyer and Pebesma 2021, to
estimate a dissimilarity index (DI) using model fitting/projection
variables, and cross-validation folds to estimate applicability limits.
See references for further information.  
The returned list includes:

- `variable_weights`, a `data.frame` of variable importance weights.

- `variable_scale_params`, a `list` of the parameters used to scale each
  variable.

- `train_DI`, a `list` containing DI estimates and thresholds per
  cross-validation type.

- `proj_aoa`, the AOA output, where 1 indicates a location within the
  AOA, 0 outside the AOA, same data type as `proj_data`.

## References

- Meyer, H., & Pebesma, E. (2021). Predicting into unknown space?
  Estimating the area of applicability of spatial prediction models.
  Methods in Ecology and Evolution, 12, 1620–1633.
  https://doi.org/10.1111/2041-210X.13650

- Meyer H, Milà C, Ludwig M, Linnenbrink J (2023). *CAST: 'caret'
  Applications for Spatial-Temporal Models*. R package,
  <https://CRAN.R-project.org/package=CAST>.

## See also

[`hib_cv_model`](https://iiasa.github.io/hibiiscus/reference/hib_cv_model.md)
[`HibCV`](https://iiasa.github.io/hibiiscus/reference/HibCV-Class.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  aoa_res <- hib_calc_aoa(mod_obj, cv_obj, proj_data)
} # } 
```
