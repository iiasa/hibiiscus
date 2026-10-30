# hib_calc_site_contrasts: Calculate differences between sites for a range of environmental variables.

This function calculates measures of differences between sites within a
study, which can then be used as covariates in models of compositional
dis/similarity. A `data.frame` (named `site_contrasts`) containing the
site-contrast pairs (`SSBS1`, `SSBS2`) and the associated 'differences'
is added to the provided `HibDatabase`.

## Usage

``` r
hib_calc_site_contrasts(
  obj,
  metrics,
  custom_metric_specs = NULL,
  baseline_lui = "Primary vegetation_Minimal use",
  overwrite = FALSE
)
```

## Arguments

- obj:

  A
  [`HibDatabase`](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
  object.

- metrics:

  A `character` vector of 'distance' metrics. Current options are;
  "lui_contrast", "geogr_dist_hav", "geogr_dist_geod" and
  "geogr_midpoint", "env_dist".

- custom_metric_specs:

  A named `list` (or `NULL`) with each named item containing the
  specifications for the named custom metric. At present, custom metrics
  must be based on `vegdist` distance metrics.  
  Essential named values within the specifications of a custom metric
  are:

  - `dist_metr` - a character indicating the metric of choice;

  - `vars` - a character vector indicating the site variables to use
    when calculating the distance metric.  

  For example:  
  `list(Gower_dist_a = list(dist_metr = "gower", vars = c("env1", "env2", "env3")), Gower_dist_b = list(dist_metr = "gower", vars = c("env4", "env5", "env6")))`  
  Metric names (as specified in the name of the sub-lists) will be
  converted to sentence-case in the output data.

- baseline_lui:

  A `character` (can be multiple entries) or `NULL`, indicating the
  baseline LUI/s against which differences are calculated per study.
  Value/s provided should match the value/s used when calculating
  compositional dis/similarity
  ([`hib_calc_div_between`](https://iiasa.github.io/hibiiscus/reference/hib_calc_div_between.md)).

- overwrite:

  A `logical` indicating whether or not any current site contrast data
  should be over-written. Default value is `FALSE`.

## Value

HibDatabase

## Details

The function currently returns `site1-site2` contrasts for all
within-study contrasts where `site1` has the `baseline_lui`. Subsequent,
filtering to ensure contrasts are based on the same effort level is
recommended and may be incrporated in future releases.

## References

https://adrianadepalma.github.io/BII_tutorial

## Examples

``` r
if (FALSE) { # \dontrun{
  obj <- hib_calc_site_contrasts(obj, 
                                 c("lui_contrast", "geogr_dist_hav"),
                                 baseline_lui = "Primary vegetation_Minimal use")
} # } 
```
