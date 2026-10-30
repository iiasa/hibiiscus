# hib_save: Save Hib\_\_\_ objects.

Saves a hibiiscus data object to disk.

## Usage

``` r
hib_save(obj, path, overwrite = FALSE)
```

## Arguments

- obj:

  A hibiiscus object, one of: `HibDatabase`, `HibLULegend`,
  `HibModelLME4`, `HibModelINLA`, `HibModelBRMS`, `HibCV`,
  `HibProjection`.

- path:

  A `character` pointing to a file (file names should have a ".rds"
  extension. Non-rds extensions will be amended to rds).

- overwrite:

  A `logical` indicating whether a previously saved file should be
  overwritten. (Default: `FALSE`)

## Details

The function saves the hibiiscus object in the specified file location.
If the object (`HibProjection`) includes `SpatRaster` data, this will be
separately saved in '.tif' file/s, with the same base-path.

## See also

[`hib_load()`](https://iiasa.github.io/hibiiscus/reference/hib_load.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  hib_save(model_obj, "an/example/file_path/hib_model.rds", T)
} # }
```
