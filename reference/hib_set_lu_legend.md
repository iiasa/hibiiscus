# hib_set_lu_legend: Set a `HibLULegend` as the LU legend in a `HibDatabase`

This function adds the `HibLULegend` to the `HibDatabase`, modifying
both the `sites` and `species` data.frames.

## Usage

``` r
hib_set_lu_legend(obj, leg)
```

## Arguments

- obj:

  A `HibDatabase`.

- leg:

  A `HibLULegend`.

## Value

`HibDatabase`

## Details

This function updates (and therefore overwrites) the previous
`Predominant_land_use`, `Use_intensity` and `LUI` columns in both the
`sites` and `species` data.frames.

## See also

[HibDatabase](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
[HibLULegend](https://iiasa.github.io/hibiiscus/reference/HibLULegend-Class.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  obj <- hib_set_lu_legend(obj, leg)
} # }
```
