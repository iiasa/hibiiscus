# Calculate between-site (within study) diversity metrics

This function calculates a range of compositional dissimilarity/distance
measures between sites (within studies) based on taxonomic diversity and
for the associated HibDatabase.

## Usage

``` r
hib_calc_div_between(
  obj,
  metrics,
  binary = FALSE,
  rel_abund = TRUE,
  rescale = TRUE,
  baseline_lui = "Primary vegetation_Minimal use"
)
```

## Arguments

- obj:

  A
  [`HibDatabase`](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
  object.

- metrics:

  Character vector of dissimilarity/distance metrics. Metrics available
  via
  [`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html):
  "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
  "jaccard", "gower", "altGower", "morisita", "horn", "mountford",
  "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord",
  "hellinger", "aitchison", "robust.aitchison". Other available metrics:
  "hypervolume", "cosine", "theta", "baselga_bray_balanced",
  "baselga_bray_gradient", "baselga_ruzicka_balanced",
  "baselga_ruzicka_gradient", "horn_dissim", "normalised_euclidean".
  Metrics can be provided in upper or lower case - internal processing
  is done using lower-case - and as as unambiguous abbreviations. Note:
  "horn" refers to the "Morisita-Horn" index in the vegan package.
  "horn_dissim" refers to (1 - Horn overlap), calculated endogenously
  within this package.

- binary:

  Logical vector - should metric be calculated using occurence values
  (`TRUE`) or abundance (`FALSE`, default). Should be a single entry, or
  of the same length as 'metrics'. Vectors of other lengths will be
  recycled to the length of 'metrics'.

- rel_abund:

  Logical vector - determining whether relative abundance should be used
  to calculate the dissimilarity metric, only applicable for some
  metrics, e.g., euclidean, hypervolume. Note: many metrics include
  rescaling of abundances as part of the calculation. Default values is
  `TRUE`.

- rescale:

  Logical vector - determining if dissimilarity measure should be
  rescaled relative to its maximum value - applies to metrics that are
  not on a 0-1 scale by default - hypervolume, euclidean (if relative
  abundance), chord, hellinger. Default values is `TRUE`.

- baseline_lui:

  Character or NULL, indicating the baseline LUI/s against which
  compositional dissimilarity is calculated per study. If NULL,
  compositional dissimilarity is calculated between all sites within
  each study. Default values is 'Primary vegetation_Minimal use' - the
  most pristine LUI within the PREDICTS database.

## Value

HibDatabase with updated '\$diversity_between' tibble

## Details

As is standard for R6 Classes, the HibDatabase is modified in-place. If
effort correction (`hib_effort_correction()`) has previously been
conducted, diversity metrics will be estimated using effort corrected
measurements.

## Examples

``` r
if (FALSE) { # \dontrun{
  obj <- hib_calc_div_between(obj, 
                              c("bray", "jaccard"),
                              c(F, T),
                              c(T, F),
                              c(F, F),
                              "Primary vegetation_Minimal use")
} # } 
```
