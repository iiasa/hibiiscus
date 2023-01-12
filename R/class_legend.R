#' Template objects for land-use legends
#' 
#' @description 
#' This is a dummy template object for the land-use legend to be used together with PREDICTS or PREDICTS-like data.
#' The idea here is to provide an object that contains a categorization of all sites in a [`Database`] object. The reason to encode
#' this information in an object is to provide a wrapper so as to load and modify existing legends (such as the default one in PREDICTS). 
#' 
#' @keywords internal, land-use
#' @import R6
#' @noRd
Legend <- R6Class("Legend",
                    # Public functions
                    public = list(
                      version = NULL,
                      last_update = NULL,
                      sites = NULL,
                      species = NULL,
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