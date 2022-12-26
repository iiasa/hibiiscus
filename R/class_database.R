#' Template objects for databases 
#' 
#' @description 
#' This is a dummy template object for a database containing PREDICTS data.
#' It contains all the functions inherent to the database. Any functions operating with the data
#' contained within the database need to be coded in here or provided separately.
#' @keywords internal
#' @import R6
#' @noRd
Database <- R6Class("Database",
                  # Public functions
                  public = list(
                    version = NULL,
                    last_update = NULL,
                    sites = NULL,
                    species = NULL,
                    initialize = function(name = NA) {
                      self$name <- name
                    },
                    print = function() {
                      cat(paste0("Hello, my name is ", self$name, ".\n"))
                    }
                  ),
                  # Private entries to be updated
                  private = list(
                    siteids = list(),
                    length = function() base::length(private$siteids)
                  )
)