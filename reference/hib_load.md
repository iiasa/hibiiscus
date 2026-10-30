# hib_load: Load Hib\_\_\_ objects.

Load a hibiiscus object from a saved RDS.

## Usage

``` r
hib_load(path)
```

## Arguments

- path:

  A `character` pointing to a file (file names should have a ".rds"
  extension, if they don't, the extension will be added).

## Value

A hibiiscus object, one of: `HibDatabase`, `HibLULegend`,
`HibModelLME4`, `HibModelINLA`, `HibModelBRMS`, `HibCV`, `HibProjection`

## Details

The function loads the hibiiscus object from the specified file
location. If the original object (HibProjection) included `SpatRaster`
data, this will be loaded from the appropriate tif file/s.

## See also

[`hib_save()`](https://iiasa.github.io/hibiiscus/reference/hib_save.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  hib_obj <- hib_load("an/example/file_path/hib_model.rds")
} # }
```
