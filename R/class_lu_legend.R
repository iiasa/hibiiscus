# Yes, I am thinking of a class for example legend_lum that has internally a 
# wrapper for recording the overall classes in terms of land-use (Primary_Forest, 
# Secondary, ...) and management (low, medium, high). The key bit is to be able 
# to link the sites (SSBS) with the entry of the legend.

# So the package should have functions to
# 
# create a new legend and save it somewhere.
# load an existing legend (PREDICTS) or a custom one.
# The object could possibly also have functions to plot the distribution of LU 
# per site / histograms.
# 
# The actual process of creating such a legend has to be done externally of course 
# and the custom legends for example for LAMASUS do not need to be part of the 
# package (in fact i would propose to store them externally).


#' Template object for land-use legend 
#' @title LU Legend object
#' @description 
#' This is the template object for a land-use legend that can be linked to a database
#' object - contains info on unique land-use and management/intensity types
#' @keywords internal
#' @import R6
#' @noRd
HibLULegend <- R6Class("HibLULegend",
                    # "lu_legend"
                    # TODO: @export
                    # Public data objects and functions
                    public = list(
                      #' @field last_update (`character()`)\cr
                      #' Identifies the most recent update to the legend object (YYYY-MM-DD)
                      last_update = NULL,
                      
                      #' @field input_data (`character()`)\cr
                      #' Identifies the source of the input data, either "file (\[...\])",
                      #' or "dataframe \[...\]".
                      input_data = NULL,
                      
                      #' @field lu_df (`data.frame`)\cr
                      #' Contains the LU legend dataframe, i.e, A data.frame 
                      #' containing columns for "SSBS", "Predominant_land_use", 
                      #' "Use_intensity" and "LUI".
                      lu_df = NULL,
                      
                      #' @field land_use (`character()`)\cr
                      #' The unique values of lu_df$Predominant_land_use.
                      land_use = NULL,
                      
                      #' @field management (`character()`)\cr
                      #' The unique values of lu_df$Use_intensity
                      management = NULL,
                      
                      #' @field land_use_management (`character()`)\cr
                      #' The unique values of lu_df$LUI
                      land_use_management = NULL,
                      
                      # Executed when created
                      #' @description
                      #' Creates a new instance of this [R6][R6::R6Class] class.
                      #' Note, this object is constructed via the [create_lu_legend]
                      #' function found within this package.
                      #' @param input (`character()` or `data.frame`)\cr
                      #'   The source data used to generate the LU legend
                      #' @param input_chr (`character()`)\cr
                      #'   Text indicating the source of the data.
                      initialize = function(input, input_chr){
                        # Check that input is either a character filepath
                        # Or dataframe/tibble
                        assertthat::assert_that(
                          (is.character(input) || is.data.frame(input)),
                          is.character(input_chr)
                          )
                        
                        if (is.character(input)){
                          assertthat::assert_that(file.exists(input))
                        }
                        
                        # ------ #
                        # Load and prepare data
                        if (is.character(input)){
                          message(crayon::green(paste0("Initializing land-use legend from file: ",
                                                       input)) )
                          lu_df <- readr::read_rds(input)
                          self$input_data <- paste0("file ('", 
                                                    input,
                                                    "')")
                            
                        }
                        if (is.data.frame(input)){
                          message(crayon::green(paste0("Initializing land-use legend from dataframe")) )
                          # : ",
                          # deparse(substitute(input))
                          
                          lu_df <- input
                          self$input_data <- paste0("dataframe ('", 
                                                    input_chr,
                                                    # deparse(substitute(input)),
                                                    "')")
                        }
                        
                        # To make subsequent data manipulation easier - e.g site_contrasts 
                        # etc - particular column names are enforced:
                        # Predominant_land_use, Use_intensity, LUI - as in PREDICTS
                        
                        # Check column names 
                        assertthat::assert_that(
                          assertthat::has_name(lu_df, "SSBS"),
                          
                          assertthat::has_name(lu_df, "Predominant_land_use"),
                          assertthat::has_name(lu_df, "Use_intensity"),
                          assertthat::has_name(lu_df, "LUI")
                        )
                        
                        # Get unique land-uses and intensities
                        un_lu <- unique(lu_df$Predominant_land_use)
                        un_i <- unique(lu_df$Use_intensity)
                        un_lui <- unique(lu_df$LUI)
                        
                        
                        # Store legend object (list of data)
                        self$last_update <- lubridate::today()
                        
                        self$lu_df <- lu_df
                        
                        self$land_use <- un_lu
                        self$management <- un_i
                        self$land_use_management <- un_lui
                        
                        # And check 
                        # Outsourced security checks
                        self$checkself()
                        
                      },
                      
                      #' @description
                      #' Self check function to ensure that the legend data.frame
                      #' contains the correct columns.
                      #'
                      #' @param df (`data.frame`).
                      #' @return (`NULL`).
                      checkself = function(df = self$lu_df){
                        # ------ #
                        # Security checks - recheck

                        assertthat::assert_that(
                          assertthat::has_name(df, "SSBS"),
                          
                          assertthat::has_name(df, "Predominant_land_use"),
                          assertthat::has_name(df, "Use_intensity"),
                          assertthat::has_name(df, "LUI")
                          
                        )
                        
                        assertthat::assert_that(
                          nrow(df) > 0,
                          length(unique(df$SSBS)) > 0
                        )
                        invisible()
                      },
                      
                      #' @description
                      #' Base printout function presents core information about the 
                      #' object.
                      #'
                      #' @return (`NULL`).
                      print = function() {
                        # Compute some stats
                        stats <- list(number_sites = length( unique(self$lu_df$SSBS) ),
                                      number_lu = length( self$land_use ),
                                      number_i = length( self$management ),
                                      number_lui = length( self$land_use_management )
                        )  
                        
                        message(
                          paste0("Hibiiscus land-use legend \n",
                                 "(Last updated: ", crayon::green( self$last_update ), ")",
                                 "\n",
                                 "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                                 "\t N land use types:              ", crayon::green(stats$number_lu), "\n",
                                 "\t N management types:            ", crayon::green(stats$number_i), "\n",
                                 "\t N land use x management types: ", crayon::green(stats$number_lui), "\n",
                                 "\n"
                          )
                        )
                        
                      },
                      
                      #' @description
                      #' Function to save the object to the specified  file-path 
                      #' ('output_path', must be '.rds').
                      #'
                      #' @param output_path (`character()`)
                      #' @return (`NULL`).
                      save_legend = function(output_path = NULL){
                        assertthat::assert_that(
                          is.character(output_path)
                        )
                        
                        if(!assertthat::has_extension(output_path, "rds")){
                          output_path <- paste0(output_path, ".rds")
                        }
                        
                        if (!dir.exists(dirname(output_path))){
                          message(crayon::yellow(paste0("\nWarning: '", dirname(output_path), "' doesn't exist.",
                                                        "\n",
                                                        "This folder will be created.\n")))
                          dir.create(dirname(output_path), recursive = T)
                        }
                        
                        if (file.exists(output_path)){
                          message(crayon::yellow(paste0("\nWarning: '", output_path, "' already exists.",
                                                        "\n",
                                                        "This file will be updated.\n")))
                        }
                        
                        message(crayon::green(paste0("Saving land-use legend to: '", output_path, "'\n")) )
                        
                        saveRDS(self, output_path)
                        assertthat::assert_that(file.exists(output_path),
                                                msg = "Writing a copy of the legend failed...")
                        invisible()
                      },
                      
                      #' @description
                      #' Function to generate figures summarising the data present
                      #' in the LU legend
                      #' 
                      #' @param x_var (`character()`) Options are: "land-use", 
                      #' "use-intensity", "land-use:use-intensity".
                      #' @return (`ggplot`).
                      plot = function(x_var = "land-use"){
                        
                        x_var <- match.arg(x_var,
                                           c("land-use", "use-intensity",
                                             "land-use:use-intensity"))
                        # Check var is in colnames(self$lu_df)
                        
                        basic_thm <- ggplot2::theme(axis.text = ggplot2::element_text(size = 16),
                                           axis.title = ggplot2::element_text(size = 20),
                                           plot.title = ggplot2::element_text(size = 20),
                                           plot.subtitle = ggplot2::element_text(size = 18),
                                           strip.text.x = ggplot2::element_text(size = 18),
                                           strip.text.y = ggplot2::element_text(size = 18),
                                           strip.background = ggplot2::element_blank(),
                                           legend.text = ggplot2::element_text(size = 16),
                                           legend.title = ggplot2::element_text(size = 18),
                                           legend.text.align = 0) 
                        
                        # Groupy by var
                        # N Sites (SSBS)
                        if (x_var == "land-use"){
                          
                          df <- self$lu_df |>
                            dplyr::group_by(Predominant_land_use) |>
                            dplyr::summarise(N = n_distinct(SSBS))
                          
                          plt <- df |>
                            ggplot2::ggplot() +
                            ggplot2::geom_hline(yintercept = 0, colour = "darkgrey", lty = 1) +
                            # geom_hline(yintercept = 100, colour = "darkgrey", lty = 2) +
                            ggplot2::geom_bar(ggplot2::aes(x = Predominant_land_use, y = N),
                                     stat = "identity",
                                     fill = "#009E73",
                                     alpha = 0.7) +
                            # geom_point(ggplot2::aes(x = Predominant_land_use, y = N),
                            #            colour = "#009E73",
                            #            alpha = 0.7,
                            #            size = 2.5) +
                            ggplot2::labs(x = "Land use",
                                 y = "Count") +
                            ggplot2::theme_bw() +
                            basic_thm + 
                            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 50,
                                                             hjust = 1,
                                                             vjust = 1),
                                  aspect.ratio = 0.75)
                          
                        }
                        
                        if (x_var == "use-intensity"){
                          
                          df <- self$lu_df |>
                            dplyr::group_by(Use_intensity) |>
                            dplyr::summarise(N = n_distinct(SSBS))
                          
                          plt <- df |>
                            ggplot2::ggplot() +
                            ggplot2::geom_hline(yintercept = 0, colour = "darkgrey", lty = 1) +
                            # geom_hline(yintercept = 100, colour = "darkgrey", lty = 2) +
                            ggplot2::geom_bar(ggplot2::aes(x = Use_intensity, y = N),
                                     stat = "identity",
                                     fill = "#009E73",
                                     alpha = 0.7) +
                            # geom_point(ggplot2::aes(x = Use_intensity, y = N),
                            #            colour = "#009E73",
                            #            alpha = 0.7,
                            #            size = 2.5) +
                            ggplot2::labs(x = "Use intensity",
                                 y = "Count") +
                            ggplot2::theme_bw() +
                            basic_thm + 
                            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 50,
                                                             hjust = 1,
                                                             vjust = 1),
                                  aspect.ratio = 0.75)
                        }
                        
                        
                        if (x_var == "land-use:use-intensity"){
                          
                          df <- self$lu_df |>
                              dplyr::group_by(Predominant_land_use, Use_intensity) |>
                              dplyr::summarise(N = n_distinct(SSBS))
                            
                          plt <- df |>
                            ggplot2::ggplot() +
                            ggplot2::geom_hline(yintercept = 0, colour = "darkgrey", lty = 1) +
                            # geom_hline(yintercept = 100, colour = "darkgrey", lty = 2) +
                            ggplot2::geom_bar(ggplot2::aes(x = Predominant_land_use, y = N,
                                         group = Use_intensity,
                                         fill = Use_intensity),
                                     position = ggplot2::position_dodge(),
                                     stat = "identity",
                                     # fill = "#009E73",
                                     alpha = 0.7) +
                            # geom_point(ggplot2::aes(x = Predominant_land_use, y = N,
                            #                group = Use_intensity,
                            #                shape = Use_intensity),
                            #            position = position_dodge(0.4),
                            #            colour = "#009E73",
                            #            alpha = 0.7,
                            #            size = 2.5) +
                            ggplot2::labs(x = "Land-use",
                                 y = "Count",
                                 # shape = "Use intensity"
                                 fill = "Use intensity") +
                            ggplot2::theme_bw() +
                            basic_thm + 
                            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 50,
                                                             hjust = 1,
                                                             vjust = 1),
                                  aspect.ratio = 0.75)
                          
                        }
                        
                        (plt)
                      }
                      
                    ),
                    # Private entries to be created
                    private = list(
                      siteids = list(),
                      length = list()
                    )
)