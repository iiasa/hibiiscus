# hib_join_sites: Join external site information to a HibDatabase object

This function joins external information to a
[`HibDatabase`](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
object. Such information could be for example additional site-specific
data that are not already present in the database.

## Usage

``` r
hib_join_sites(obj, site_info)
```

## Arguments

- obj:

  A `HibDatabase` object.

- site_info:

  A `data.frame` or similar of categorical or numeric information to be
  added to the site dataframe within the `HibDatabase`. Must share at
  least one column with `HibDatabase$sites`, e.g., 'SSBS'.

## Value

HibDatabase

## Examples

``` r
if (FALSE) { # \dontrun{
  hib_db <- hib_join_sites(hib_db, new_site_data_df)
} # } 
```
