# prep_globiom: Prepare GLOBIOM output for projection.

Function to format a prepared GLOBIOM netCDF file for use in Hibiiscus
projections.

## Usage

``` r
prep_globiom(
  file_path,
  out_as_rast = TRUE,
  incl_pixel_area = TRUE,
  vars_to_include = NULL,
  vars_to_ignore = NULL
)
```

## Arguments

- file_path:

  A filepath (`character`) pointing to a GLOBIOM output in netCDF
  format.

- out_as_rast:

  A `logical` denoting if the output type should be a (list) of
  `SpatRaster`, or a `stars` object (default: `TRUE`, i.e. returns
  `SpatRaster`).

- incl_pixel_area:

  A [`logical`](https://rdrr.io/r/base/logical.html) denoting if the
  output should include the 'pixel_area_ha' variable (Default: `TRUE`).

- vars_to_include:

  A `character` vector of variables to be extracted. If `NULL`, all
  variables will be extracted.

- vars_to_ignore:

  A `character` vector of variables to be ignored If `NULL`, all
  variables will be extracted - unless vars_to_include is not `NULL`.

## Value

A (list of) `SpatRaster` stacks with the formatted GLOBIOM predictors.

## Details

This function expects a downscaled GLOBIOM output. Variables are
extracted and converted to raster format, with each raster stack
representing a particular year. Likely of little use for anyone outside
IIASA.

## Examples

``` r
if (FALSE) { # \dontrun{
# Expects a filename pointing to a globiom-netCDF file.
globiom_rast <- prep_globiom("path/to/globiom_data.nc")
} # }
```
