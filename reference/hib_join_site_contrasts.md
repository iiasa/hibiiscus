# hib_join_site_contrasts: Join external site-contrast information to a HibDatabase object

This function joins external information to a
[`HibDatabase`](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
object. Such information could be for example additional
site-contrast-specific data that are not already present in the
database.

## Usage

``` r
hib_join_site_contrasts(obj, site_contrast_info)
```

## Arguments

- obj:

  A `HibDatabase` object.

- site_contrast_info:

  A `data.frame` or similar of categorical or numeric information to be
  added to the database. Must share at least one column with
  `HibDatabase$site_contrasts`. Ideally, both/one of 'SSBS1' and 'SSBS2'
  are present.

## Value

HibDatabase

## Examples

``` r
if (FALSE) { # \dontrun{
  hib_db <- hib_join_site_contrasts(hib_db, new_site_contrast_data_df)
} # } 
```
