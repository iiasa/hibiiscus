#' Template objects for databases 
#' 
#' @description 
#' This is the template object for a database containing PREDICTS or PREDICTS-like data.
#' It contains all the functions inherent to the database. Any functions in the package that is 
#' operating or interacting with the data (preparation, linking, modelling) need to be coded in here 
#' or provided separately.
#' @keywords internal
#' @import R6
#' @noRd
Database <- R6Class("Database",
                  # Public data objects and functions
                  public = list(
                    version = NULL,
                    last_update = NULL,
                    paths = NULL,
                    sites = NULL,
                    species = NULL,
                    # Executed when created
                    initialize = function(version, paths){
                      assertthat::assert_that(is.character(version), is.list(paths))

                      # ------ #
                      message(crayon::green(paste0("Initializing input datasets for version ", version, "...")) )
                      
                      # Now load both datasets
                      if(version == "2016"){
                        sites <- readr::read_rds(paths$sites)
                        species <- readr::read_rds(paths$species)
                      }
                      # Outsourced security checks
                      self$checkself(sites, species)
                      
                      # Save inputs in list
                      self$version <- version
                      self$last_update <- lubridate::today()
                      self$paths <- paths
                      self$sites <- sites
                      self$species <- species
                    },
                    # Self check function that everything is in order
                    checkself = function(sites = self$sites, species = self$species){
                      # ------ #
                      # Security checks
                      assertthat::assert_that(
                        assertthat::has_name(sites, "Source_ID"),
                        assertthat::has_name(sites, "SS"),
                        assertthat::has_name(sites, "Diversity_metric"),
                        assertthat::has_name(sites, "SSBS"),
                        assertthat::has_name(sites, "Predominant_land_use"),
                        # And the species data
                        assertthat::has_name(species, "Source_ID"),
                        assertthat::has_name(species, "SS"),
                        assertthat::has_name(species, "SSBS"),
                        assertthat::has_name(species, "Best_guess_binomial"),
                        assertthat::has_name(species, "Measurement")
                      )
                      # Finally, check that all sites are in species
                      assertthat::assert_that(
                        nrow(sites)>0, nrow(species)>0,
                        all( unique(sites$Source_ID) %in% unique(species$Source_ID) )
                      )
                      invisible()
                    },
                    # Base printout function
                    print = function() {
                      # Compute some stats
                      stats <- list(number_sites = nrow(self$sites),
                                    number_lu = length( unique(self$sites$Predominant_land_use) ),
                                    number_lui = length(unique(self$sites$Use_intensity)),
                                    number_species = length( unique(self$species$Best_guess_binomial) ))
                      message(
                        paste0("Hibiiscus PREDICTS database (including versions: ", crayon::bold(self$version)," )","\n",
                               "(Last updated: ", crayon::green( self$last_update ), ")",
                               "\n",
                               "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                               "\t Land use (class x intensity):  ", crayon::green(stats$number_lu), " x ",
                                                                     crayon::green(stats$number_lui), "\n",
                               "\t Number of species:             ", crayon::green(stats$number_species)
                               )
                        )
                    }
                  ),
                  # Private entries to be created
                  private = list(
                    siteids = list(),
                    length = function() base::length(self$sites$SS)
                  )
)