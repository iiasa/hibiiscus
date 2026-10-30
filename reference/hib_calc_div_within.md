# Calculate within-site diversity metrics

This function calculates multiple, user specified measures of site-level
(taxonomic) diversity for sites in the associated HibDatabase.

## Usage

``` r
hib_calc_div_within(obj, metrics)
```

## Arguments

- obj:

  A
  [`HibDatabase`](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
  object.

- metrics:

  Character vector, any/all of; "Total_abundance", "Species_richness",
  "Chao1_diversity", "Simpson_diversity", "Rarefied_species_richness",
  "Mean_abundance", "Geometric_mean_abundance", "Shannon_diversity",
  "HillShannon_diversity", "HillSimpson_diversity", "PIE", "Equalness".
  Metrics can be provided in upper or lower case and as unambiguous
  abbreviations (internal processing is done using lower-case).

## Value

HibDatabase with updated '\$diversity_within' tibble

## Details

As is standard for R6 Classes, the HibDatabase is modified in-place. If
effort correction (`hib_effort_correction()`) has previously been
conducted, diversity metrics will be estimated using effort corrected
measurements.

## Examples

``` r
if (FALSE) { # \dontrun{
  hib_db <- hib_calc_div_within(hib_db, c("Total_abundance", "Species_richness"))
} # }
 
```
