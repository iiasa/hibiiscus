# hib_update_db: Update an existing `HibDatabase` with further data

This function updates an existing `HibDatabase`, for example by adding
the 2022 version of the PREDICTS data to the 2016 version.

## Usage

``` r
hib_update_db(obj, path_species, path_sites = "", dataset_name = "")
```

## Arguments

- obj:

  A
  [`HibDatabase`](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
  object.

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

Newly provided data (from files) will be appended to the previously
created `HibDatabase`. Species data will be added to the species
data.frame, and site data to the site data.frame. Meta-data/descriptors
in the `HibDatabase` will be updated to reflect the additions. Note:
Prior filtering/processing of the initial `HibDatabase` will not be
applied to the newly provided data. Please load all data of interest
into the `HibDatabase` prior to further processing.

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
  hibDB <- hib_create_db(path_species, path_sites, "predicts_2016")
  path_species_upd <- "rawdata/PREDICTS/2022/database2ndrelease.rds"
  hibDB <- update_hib(hibDB, path_species_upd, "", "predicts_2022")
  hibDB
} # } 
```
