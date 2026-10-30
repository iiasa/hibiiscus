# hib_create_lu_legend: Create a 'new' land-use legend.

This function creates a
[`R6`](https://r6.r-lib.org/reference/R6Class.html)
[`HibLULegend`](https://iiasa.github.io/hibiiscus/reference/HibLULegend-Class.md)
object containing the provided land-use legend (i.e., linking sites,
SSBS, to land-uses/intensities).

## Usage

``` r
hib_create_lu_legend(input)
```

## Arguments

- input:

  Either a `character` pointing to a saved rds file, or a
  `data.frame`-like object within the environment.

## Value

`HibLULegend`

## Details

The provided data.frame must include the following columns: `SSBS`,
`Predominant_land_use`, `Use_intensity` and `LUI`.

## See also

[`HibLULegend`](https://iiasa.github.io/hibiiscus/reference/HibLULegend-Class.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  lu_leg <- hib_create_lu_legend("../Data/landuse_dataframe.rds")
} # }
```
