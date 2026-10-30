# HibLULegend

This is the template object for an
[R6](https://r6.r-lib.org/reference/R6Class.html) land-use legend that
can be linked to a
[HibDatabase](https://iiasa.github.io/hibiiscus/reference/HibDatabase-Class.md)
object. The legend should contain information on land-use and
management/intensity types for sites with the database.

## Details

Template object for land-use legend

## Public fields

- `last_update`:

  ([`character()`](https://rdrr.io/r/base/character.html))  
  Identifies the most recent update to the legend object (YYYY-MM-DD)

- `input_data`:

  ([`character()`](https://rdrr.io/r/base/character.html))  
  Identifies the source of the input data, either "file (\[...\])", or
  "dataframe \[...\]".

- `lu_df`:

  (`data.frame`)  
  Contains the LU legend dataframe, i.e, A data.frame containing columns
  for "SSBS", "Predominant_land_use", "Use_intensity" and "LUI".

- `land_use`:

  ([`character()`](https://rdrr.io/r/base/character.html))  
  The unique values of lu_df\$Predominant_land_use.

- `management`:

  ([`character()`](https://rdrr.io/r/base/character.html))  
  The unique values of lu_df\$Use_intensity

- `land_use_management`:

  ([`character()`](https://rdrr.io/r/base/character.html))  
  The unique values of lu_df\$LUI

## Methods

### Public methods

- [`HibLULegend$new()`](#method-HibLULegend-initialize)

- [`HibLULegend$checkself()`](#method-HibLULegend-checkself)

- [`HibLULegend$print()`](#method-HibLULegend-print)

- [`HibLULegend$plot()`](#method-HibLULegend-plot)

- [`HibLULegend$clone()`](#method-HibLULegend-clone)

------------------------------------------------------------------------

### `HibLULegend$new()`

Creates a new instance of this
[R6](https://r6.r-lib.org/reference/R6Class.html) class. Note, this
object is constructed via the
[hib_create_lu_legend](https://iiasa.github.io/hibiiscus/reference/hib_create_lu_legend.md)
function found within this package.

#### Usage

    HibLULegend$new(input, input_chr)

#### Arguments

- `input`:

  ([`character()`](https://rdrr.io/r/base/character.html) (rds filepath)
  or `data.frame`)  
  The source data used to generate the LU legend

- `input_chr`:

  ([`character()`](https://rdrr.io/r/base/character.html))  
  Text indicating the source of the data.

#### Returns

(`HibLULegend`)  

------------------------------------------------------------------------

### `HibLULegend$checkself()`

Self check function to ensure that the legend data.frame contains the
correct columns.

#### Usage

    HibLULegend$checkself(df = self$lu_df)

#### Arguments

- `df`:

  (`data.frame`).

#### Returns

(`NULL`).

------------------------------------------------------------------------

### `HibLULegend$print()`

Base printout function presents core information about the object.

#### Usage

    HibLULegend$print()

#### Returns

(`NULL`).

------------------------------------------------------------------------

### `HibLULegend$plot()`

Function to generate figures summarising the data present in the LU
legend

#### Usage

    HibLULegend$plot(x_var = "land-use")

#### Arguments

- `x_var`:

  ([`character()`](https://rdrr.io/r/base/character.html)) Options are:
  "land-use", "use-intensity", "land-use:use-intensity".

#### Returns

(`ggplot`).

------------------------------------------------------------------------

### `HibLULegend$clone()`

The objects of this class are cloneable with this method.

#### Usage

    HibLULegend$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
