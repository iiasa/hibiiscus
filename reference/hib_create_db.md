# hib_create_db: Create a `HibDatabase` object from input files

This function creates an
[`R6`](https://r6.r-lib.org/reference/R6Class.html) `HibDatabase` object
containing the loaded `PREDICTS`-like database.

## Usage

``` r
hib_create_db(path_species, path_sites = "", dataset_name = "")
```

## Arguments

- path_species:

  A character pointing to a `PREDICTS` formatted abundance/ occurrence
  database.

- path_sites:

  Optional character argument, pointing to a `PREDICTS` formatted site
  database. Default value is `""`. If `""` is provided, site data will
  be generated from the species data.

- dataset_name:

  A `character` indicating the dataset used to create the HibDatabase.
  Default is `""`. If `""` is supplied, the parameter will be internally
  updated to the filename associated with `path_species`.

## Value

An object of class `HibDatabase`, containing site and species
information which can be further modified for subsequent ecological
analysis.

## Details

This function assumes that a species database file is explicitly
provided, with a site database file being optional. Provided file paths
are expected to include a '.rds' extension.  
All files are loaded and compiled into a single database object that can
be stored for future use (see
[`hib_save()`](https://iiasa.github.io/hibiiscus/reference/hib_save.md)).
This object is then used in other functions of this package to prepare
or filter the database and enable ecological impact analyses.

## References

- 1.  Hudson, L.N., Newbold, T., Contu, S., Hill, S.L.L., Lysenko, I.,
      De Palma, A., Phillips, H.R.P., Senior, R.A., Bennett, D.J.,
      Booth, H., et al. (2014). The PREDICTS database: a global database
      of how local terrestrial biodiversity responds to human impacts.
      Ecol. Evol. 4, 4701–4735. 10.1002/ece3.1303.

- 1.  Hudson, L.N., Newbold, T., Contu, S., Hill, S.L.L., Lysenko, I.,
      De Palma, A., Phillips, H.R.P., Alhusseini, T.I., Bedford, F.E.,
      Bennett, D.J., et al. (2017). The database of the PREDICTS
      (Projecting Responses of Ecological Diversity In Changing
      Terrestrial Systems) project. Ecol. Evol. 7, 145–188.
      10.1002/ece3.2579.

- 1.  Lawrence Hudson; Tim Newbold; Sara Contu; Samantha L L Hill et al.
      (2023). The 2016 release of the PREDICTS database V1.1. Natural
      History Museum. https://doi.org/10.5519/j4sh7e0w

- 1.  Sara Contu; Adriana De Palma; Rachel Bates; Jessica Borer et al.
      (2022). Release of data added to the PREDICTS database (November
      2022). Natural History Museum. https://doi.org/10.5519/jg7i52dg

## Examples

``` r
if (FALSE) { # \dontrun{
  path_species <- "rawdata/PREDICTS/2016/database.rds"
  path_sites <- "rawdata/PREDICTS/2016/sites.rds"
  hibDB <- hib_create_db(path_species, path_sites)
  hibDB
} # }
```
