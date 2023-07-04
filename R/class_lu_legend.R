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
#' 
#' @description 
#' This is the template object for a land-use legend that can be linked to a database
#' object - contains info on unique land-use and management/intensity types
#' @keywords internal
#' @import R6
#' @noRd
lu_legend <- R6Class("lu_legend",
                    # Public data objects and functions
                    public = list(
                      # version = NULL,
                      last_update = NULL,
                      
                      input_data = NULL,
                      lu_df = NULL,
                      land_use = NULL,
                      management = NULL,
                      land_use_management = NULL,
                      # processing_flags = NULL,
                      
                      # Executed when created
                      initialize = function(input){
                        # Check that input is either a character filepath
                        # Or dataframe/tibble
                        assertthat::assert_that(
                          (is.character(input) || is.data.frame(input))
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
                          self$input_data <- input
                            
                        }
                        if (is.data.frame(input)){
                          message(crayon::green(paste0("Initializing land-use legend from dataframe")) )
                          # : ",
                          # deparse(substitute(input))
                          
                          lu_df <- input
                          self$input_data <- "dataframe" #deparse(substitute(input))
                        }
                        
                        
                        # print(lu_df)
                        
                        # TODO:
                        # To make subsequent data manipulation easier - e.g site_contrasts etc
                        # It might be better to enforce particular column names
                        # E.g., Predominant_land_use, Use_intensity, LUI - as in PREDICTS?
                        
                        # Check column names - SSBS, LU or similar, Management or similar, LUI or similar
                        assertthat::assert_that(
                          assertthat::has_name(lu_df, "SSBS"),
                          # any(c("Predominant_land_use", "Predominant_habitat", "Land_use") %in% colnames(lu_df)),
                          # any(c("Management_intensity", "Use_intensity") %in% colnames(lu_df)),
                          # any(c("Land_use_intensity", "LUI") %in% colnames(lu_df))
                          # (assertthat::has_name(lu_df, "Predominant_land_use") || assertthat::has_name(lu_df, "Predominant_habitat") || assertthat::has_name(lu_df, "Land_use")),
                          # (assertthat::has_name(lu_df, "Management_intensity") || assertthat::has_name(lu_df, "Use_intensity")),
                          # (assertthat::has_name(lu_df, "Land_use_intensity") || assertthat::has_name(lu_df, "LUI"))
              
                          assertthat::has_name(lu_df, "Predominant_land_use"),
                          assertthat::has_name(lu_df, "Use_intensity"),
                          assertthat::has_name(lu_df, "LUI")
                        )
                        
                        # Get unique land-uses and intensities
                        # lu_col <- c("Predominant_land_use", "Predominant_habitat", "Land_use")
                        # lu_col <- lu_col[lu_col %in% colnames(lu_df)]
                        # assertthat::assert_that(
                        #   length(lu_col) == 1
                        # )
                        # 
                        # i_col <- c("Management_intensity", "Use_intensity")
                        # i_col <- i_col[i_col %in% colnames(lu_df)]
                        # assertthat::assert_that(
                        #   length(i_col) == 1
                        # )
                        # 
                        # lui_col <- c("Land_use_intensity", "LUI")
                        # lui_col <- lui_col[lui_col %in% colnames(lu_df)]
                        # assertthat::assert_that(
                        #   length(lui_col) == 1
                        # )
                        
                        # un_lu <- unique(lu_df[[lu_col]])
                        # un_i <- unique(lu_df[[i_col]])
                        # un_lui <- unique(lu_df[[lui_col]])
                        
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
                        
                        # Add some data processing flags
                        # self$processing_flags <- list(measures_effort_corrected = F,
                        #                               abundance_rescaled = F,
                        #                               filtering_steps = list())
                      },
                      # Self check function that everything is in order
                      checkself = function(df = self$lu_df){
                        # ------ #
                        # Security checks - recheck

                        assertthat::assert_that(
                          assertthat::has_name(df, "SSBS"),
                          # any(c("Predominant_land_use", "Predominant_habitat", "Land_use") %in% colnames(df)),
                          # any(c("Management_intensity", "Use_intensity") %in% colnames(df)),
                          # any(c("Land_use_intensity", "LUI") %in% colnames(df))
                          # assertthat::has_name(lu_df, "Predominant_land_use") || assertthat::has_name(lu_df, "Predominant_habitat") || assertthat::has_name(lu_df, "Land_use"),
                          # assertthat::has_name(lu_df, "Management_intensity") || assertthat::has_name(lu_df, "Use_intensity"),
                          # assertthat::has_name(lu_df, "Land_use_management") || assertthat::has_name(lu_df, "LUM") ||assertthat::has_name(lu_df, "Land_use_intensity") || assertthat::has_name(lu_df, "LUI")
                          assertthat::has_name(df, "Predominant_land_use"),
                          assertthat::has_name(df, "Use_intensity"),
                          assertthat::has_name(df, "LUI")
                          
                        )
                        # Finally, check that all sites are in species
                        # NOTE RC - source ID is just the paper...
                        assertthat::assert_that(
                          nrow(df) > 0,
                          length(unique(df$SSBS)) > 0
                        )
                        invisible()
                      },
                      # Base printout function
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
                      # Save object
                      save_legend = function(output_path = NULL){
                        assertthat::assert_that(
                          is.character(output_path)
                        )
                        
                        if(!assertthat::has_extension(output_path, "rds")){
                          output_path <- paste0(output_path, ".rds")
                        }
                        
                        if (!dir.exists(dirname(output_path))){
                          message(crayon::yellow(paste0("Warning: '", dirname(output_path), "' doesn't exist.",
                                                        "\n",
                                                        "This folder will be created.\n")))
                          dir.create(dirname(output_path), recursive = T)
                        }
                        
                        if (file.exists(output_path)){
                          message(crayon::yellow(paste0("Warning: '", output_path, "' already exists.",
                                                        "\n",
                                                        "This file will be updated.\n")))
                        }
                        
                        message(crayon::green(paste0("Saving land-use legend to: '", output_path, "' ...\n")) )
                        
                        saveRDS(self, output_path)
                        assertthat::assert_that(file.exists(output_path),
                                                msg = "Writing a copy of the legend failed...")
                        invisible()
                      },
                      
                      # TODO:
                      # summary plots - distr of luis... - either lu, i, or lui
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
                          # lu_col <- c("Predominant_land_use", "Predominant_habitat", "Land_use")
                          # lu_col <- lu_col[lu_col %in% colnames(self$lu_df)]
                          # 
                          # df <- self$lu_df %>%
                          #   dplyr::group_by(!!sym(lu_col)) %>%
                          #   dplyr::summarise(N = n_distinct(SSBS))
                          # print(head(df))
                          
                          df <- self$lu_df %>%
                            dplyr::group_by(Predominant_land_use) %>%
                            dplyr::summarise(N = n_distinct(SSBS))
                          
                          plt <- df %>%
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
                          # i_col <- c("Management_intensity", "Use_intensity")
                          # i_col <- i_col[i_col %in% colnames(self$lu_df)]
                          # 
                          # df <- self$lu_df %>%
                          #   dplyr::group_by(!!sym(i_col)) %>%
                          #   dplyr::summarise(N = n_distinct(SSBS))
                          # print(head(df))
                          
                          df <- self$lu_df %>%
                            dplyr::group_by(Use_intensity) %>%
                            dplyr::summarise(N = n_distinct(SSBS))
                          
                          plt <- df %>%
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
                          
                          # lu_col <- c("Predominant_land_use", "Predominant_habitat", "Land_use")
                          # lu_col <- lu_col[lu_col %in% colnames(self$lu_df)]
                          # 
                          # i_col <- c("Management_intensity", "Use_intensity")
                          # i_col <- i_col[i_col %in% colnames(self$lu_df)]
                          # 
                          # df <- self$lu_df %>%
                          #   dplyr::group_by(!!sym(lu_col), !!sym(i_col)) %>%
                          #   dplyr::summarise(N = n_distinct(SSBS))
                          # print(head(df))
                          
                          df <- self$lu_df %>%
                              dplyr::group_by(Predominant_land_use, Use_intensity) %>%
                              dplyr::summarise(N = n_distinct(SSBS))
                            
                          plt <- df %>%
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
                        # %%%
                        # lu_col <- c("Predominant_land_use", "Predominant_habitat", "Land_use")
                        # lu_col <- lu_col[lu_col %in% colnames(lu_df)]
                        # assertthat::assert_that(
                        #   length(lu_col == 1)
                        # )
                        # 
                        # i_col <- c("Management_intensity", "Use_intensity")
                        # i_col <- i_col[i_col %in% colnames(lu_df)]
                        # assertthat::assert_that(
                        #   length(i_col == 1)
                        # )
                        # 
                        # lui_col <- c("Land_use_intensity", "LUI")
                        # lui_col <- lui_col[lui_col %in% colnames(lu_df)]
                        # assertthat::assert_that(
                        #   length(lui_col == 1)
                        # )
                        # %%%
                        
                        # plot
                        # %%%
                        # df %>%
                        #   ggplot() +
                        #   geom_hline(yintercept = 0, colour = "darkgrey", lty = 1) +
                        #   geom_hline(yintercept = 100, colour = "darkgrey", lty = 2) +
                        #   geom_point(aes(x = var, y = N,
                        #                  group = Use_intensity,
                        #                  shape = Use_intensity),
                        #              size = 2,
                        #              position = position_dodge(0.4)) +
                        #   geom_text(aes(x = var,
                        #                 y = -100,
                        #                 label = N_tot),
                        #             data = df %>%
                        #               group_by(var) %>%
                        #               summarise(N_tot = sum(N)) ) +
                        #   labs(x = "Land use",
                        #        y = "Count",
                        #        shape = "Use intensity") +
                        #   theme_bw() +
                        #   basic_thm + 
                        #   theme(axis.text.x = element_text(angle = 50,
                        #                                    hjust = 1,
                        #                                    vjust = 1),
                        #         aspect.ratio = 0.75)
                        # # %%%
                        
                        (plt)
                      }
                      
                    ),
                    # Private entries to be created
                    private = list(
                      siteids = list(),
                      length = list()
                    )
)