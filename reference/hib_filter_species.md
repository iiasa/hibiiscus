# hib_filter_species: Filter the database object by species-linked features.

This function returns a database object containing only species (and
associated sites) whose characteristics are retained.

## Usage

``` r
hib_filter_species(obj, condition)
```

## Arguments

- obj:

  A
  [`HibDatabase`](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
  object.

- condition:

  An expression identifying how to filter the species data.frame, as
  would be passed to
  [`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html).

## Value

`HibDatabase`

## Examples

``` r
if (FALSE) { # \dontrun{
  hib_obj <- hib_filter_species(hib_obj, Class == "Insecta")
} # }
```
