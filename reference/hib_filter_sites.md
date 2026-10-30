# hib_filter_sites: Filter the database object by site-linked variables

This function returns a database object containing only sites (and
associated species) whose characteristics are retained.

## Usage

``` r
hib_filter_sites(obj, condition)
```

## Arguments

- obj:

  A
  [`HibDatabase`](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
  object.

- condition:

  An expression identifying how to filter the sites data.frame, as would
  be passed to
  [`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html).

## Value

`HibDatabase`

## Examples

``` r
if (FALSE) { # \dontrun{
  hib_obj <- hib_filter_sites(hib_obj, UN_region == "Europe")
} # }
```
