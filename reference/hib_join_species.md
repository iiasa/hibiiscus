# hib_join_species: Join external species information to a HibDatabase object

This function joins external information to a
[`HibDatabase`](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
object. Such information could be for example prepared data on species
traits.

## Usage

``` r
hib_join_species(obj, spp_info)
```

## Arguments

- obj:

  A `HibDatabase` object.

- spp_info:

  A `data.frame` or similar of categorical or numeric information to be
  added to the species data in the `HibDatabase`. Must share at least
  one column with `HibDatabase$species`, e.g. 'Best_guess_binomial'.

## Value

HibDatabase

## Examples

``` r
if (FALSE) { # \dontrun{
  hib_db <- hib_join_species(hib_db, new_spp_data_df)
} # } 
```
