#' Template object for land-use legend 
#' @title HibLULegend
#' @description 
#' This is the template object for an [R6][R6::R6Class] land-use legend that can 
#' be linked to a [HibDatabase] object. The legend should contain information on 
#' land-use and management/intensity types for sites with the database.
#' @keywords classes
#' @import R6
#' @rdname HibLULegend-Class
#' @export
HibLULegend <- R6Class("HibLULegend",
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
                         #' Note, this object is constructed via the [hib_create_lu_legend]
                         #' function found within this package.
                         #' @param input (`character()` (rds filepath) or `data.frame`)\cr
                         #'   The source data used to generate the LU legend
                         #' @param input_chr (`character()`)\cr
                         #'   Text indicating the source of the data.
                         #'   
                         #' @return (`HibLULegend`)\cr
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
                           # HibLULegend 
                           if (is.character(input)){
                             message(crayon::green(paste0("\nInitialising HibLULegend from file: ",
                                                          input, "...\n")) )
                             lu_df <- readRDS(input)
                             self$input_data <- paste0("file ('", 
                                                       input,
                                                       "')")
                             
                           }
                           if (is.data.frame(input)){
                             message(crayon::green(paste0("\nInitialising HibLULegend from data.frame: ",
                                                          input_chr, "...\n")) )
                             
                             lu_df <- input
                             self$input_data <- paste0("data.frame ('", 
                                                       input_chr,
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
                           self$last_update <- Sys.time()
                           
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
                                    "Last updated: ", crayon::green( self$last_update ),
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
                         #' Function to generate figures summarising the data present
                         #' in the LU legend
                         #' 
                         #' @param x_var (`character()`) Options are: "land-use", 
                         #' "use-intensity", "land-use:use-intensity".
                         #' @return (`ggplot`).
                         plot = function(x_var = "land-use"){
                           rlang::check_installed("ggplot2", reason = "to plot hibiiscus objects.")
                           
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
                                                 alpha = 0.7) +
                               ggplot2::labs(x = "Land-use",
                                             y = "Count",
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