# HibDatabase

This is the template object for an
[R6](https://r6.r-lib.org/reference/R6Class.html) database containing
PREDICTS or PREDICTS-like data. Functions are provided to prepare the
data for modelling and visualise the data.

## Details

Template objects for databases

## References

De Palma, A., Abrahamczyk, S., Aizen, M. et al. Predicting bee community
responses to land-use changes: Effects of geographic and taxonomic
biases. Sci Rep 6, 31153 (2016). https://doi.org/10.1038/srep31153

## Public fields

- `dataset_name`:

  ([`character()`](https://rdrr.io/r/base/character.html))  
  Identifies the dataset used.

- `last_update`:

  ([`character()`](https://rdrr.io/r/base/character.html))  
  Identifies the most recent update to the Database (YYYY-MM-DD).

- `paths`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  Identifies the files used to create the Database object.

- `sites`:

  (`data.frame`)  
  Stores information on the sites within the Database.

- `species`:

  (`data.frame`)  
  Stores information on the species within the Database.

- `diversity_within`:

  (`tibble`)  
  Stores information on biodiversity metrics calculated for each site
  (i.e. alpha diversity).

- `diversity_between`:

  (`tibble`)  
  Stores information on biodiversity metrics calculated between sites,
  within studies (i.e. beta diversity).

- `site_contrasts`:

  (`data.frame`)  
  Stores information on 'distances' between sites, within studies (e.g.,
  geographic distance).

- `processing_flags`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  Stores information on processing that has occurred (e.g., filtering
  data, adding data).

## Methods

### Public methods

- [`HibDatabase$new()`](#method-HibDatabase-initialize)

- [`HibDatabase$checkself()`](#method-HibDatabase-checkself)

- [`HibDatabase$print()`](#method-HibDatabase-print)

- [`HibDatabase$hib_effort_correction()`](#method-HibDatabase-hib_effort_correction)

- [`HibDatabase$plot()`](#method-HibDatabase-plot)

- [`HibDatabase$clone()`](#method-HibDatabase-clone)

------------------------------------------------------------------------

### `HibDatabase$new()`

Creates a new instance of this
[R6](https://r6.r-lib.org/reference/R6Class.html) class. Note, this
object is constructed via the
[hib_create_db](https://iiasa.github.io/hibiiscus/reference/hib_create_db.md)
and
[hib_update_db](https://iiasa.github.io/hibiiscus/reference/hib_update_db.md)
functions found within this package.

#### Usage

    HibDatabase$new(dataset_name, paths)

#### Arguments

- `dataset_name`:

  ([`character()`](https://rdrr.io/r/base/character.html))  
  Text indicating the dataset contributing to the database being
  created.

- `paths`:

  ([`list()`](https://rdrr.io/r/base/list.html))  
  List of paths to the 'species' and 'site' data to be loaded.

#### Returns

(`HibDatabase`)  

------------------------------------------------------------------------

### `HibDatabase$checkself()`

Self check function to ensure that the site and species data.frames
contain the correct columns and align in content.

#### Usage

    HibDatabase$checkself(sites = self$sites, species = self$species)

#### Arguments

- `sites`:

  (`data.frame`).

- `species`:

  (`data.frame`).

#### Returns

(`NULL`).

------------------------------------------------------------------------

### `HibDatabase$print()`

Base printout function presents core information about the Database
object.

#### Usage

    HibDatabase$print()

#### Returns

(`NULL`).

------------------------------------------------------------------------

### `HibDatabase$hib_effort_correction()`

Function to apply effort correction to the species data. Uses a simple,
linear scaling, such that sampling effort of a site is scaled relative
to the maximum value per study, and multiplied by the original measure
of diversity.

#### Usage

    HibDatabase$hib_effort_correction()

#### Returns

(`HibDatabase`).

------------------------------------------------------------------------

### `HibDatabase$plot()`

Function to generate figures summarising the site-level data present in
the HibDatabase.

#### Usage

    HibDatabase$plot(type = "map-sites")

#### Arguments

- `type`:

  ([`character()`](https://rdrr.io/r/base/character.html)) Options are:
  "map-sites", "map-studies", "biome-N", "ecoregion-N", "realm-N",
  "class-N".

#### Returns

(`ggplot`).

------------------------------------------------------------------------

### `HibDatabase$clone()`

The objects of this class are cloneable with this method.

#### Usage

    HibDatabase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
