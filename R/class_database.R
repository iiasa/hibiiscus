#' Template objects for databases 
#' @title HibDatabase Object
#' @description 
#' This is the template object for a database containing PREDICTS or PREDICTS-like data.
#' It contains all the functions inherent to the database. Any functions in the package that is 
#' operating or interacting with the data (preparation, linking, modelling) need to be coded in here 
#' or provided separately.
#' @keywords internal
#' @import R6
#' @noRd
HibDatabase <- R6Class("HibDatabase",
                  
                  # Public data objects and functions
                  public = list(
                    #' @field version (`character()`)\cr
                    #' Identifies the version of data used.
                    version = NULL,
                    
                    #' @field last_update (`character()`)\cr
                    #' Identifies the most recent update to the Database (YYYY-MM-DD).
                    last_update = NULL,
                    
                    #' @field paths (`list()`)\cr
                    #' Identifies the files used to create the Database object.
                    paths = NULL,
                    
                    #' @field sites (`data.frame`)\cr
                    #' Stores information on the sites within the Database.
                    sites = NULL,
                    
                    #' @field species (`data.frame`)\cr
                    #' Stores information on the species within the Database.
                    species = NULL,
                    
                    #' @field diversity_within (`tibble`)\cr
                    #' Stores information on biodiversity metrics calculated for
                    #' each site (i.e. alpha diversity).
                    diversity_within = NULL,
                    
                    #' @field diversity_between (`tibble`)\cr
                    #' Stores information on biodiversity metrics calculated 
                    #' between sites, within studies (i.e. beta diversity).
                    diversity_between = NULL,
                    
                    #' @field site_contrasts (`data.frame`)\cr
                    #' Stores information on 'distances' between sites, within 
                    #' studies (e.g., geographic distance).
                    site_contrasts = NULL,
                    
                    #' @field processing_flags (`list()`)\cr
                    #' Stores information on processing that has occurred (e.g.,
                    #' filtering data, adding data).
                    processing_flags = NULL,


                    # Executed when created
                    #' @description
                    #' Creates a new instance of this [R6][R6::R6Class] class.
                    #' Note, this object is constructed via the [create_hib] and 
                    #' [update_hib] functions found within this package.
                    #' @param version (`character()`)\cr
                    #'   Text indicating the version of the database being created.
                    #' @param paths (`list()`)\cr
                    #'   List of paths to the 'species' and 'site' data to be loaded.
                    initialize = function(version, paths){
                      assertthat::assert_that(
                        is.character(version), 
                        is.list(paths)
                        )
                      
                      # path$sites/species is character - already checked in create_hib
                      # ------ #

                      message(crayon::green(paste0("Initializing input datasets for ", 
                                                   version, " data ...\n")) )
                      
                      # %%%%%
                      # Load species data
                      assertthat::assert_that(
                        assertthat::has_extension(paths$species, "rds")
                        )

                      species <- readr::read_rds(paths$species)

                      species$version <- version
                      
                      # If predicts2022 column names, rename to align with 2016...
                      if ("Predominant_habitat" %in% colnames(species)){
                        species <- species |>
                          dplyr::rename(Predominant_land_use = Predominant_habitat)
                      }
                      if ("Source_for_predominant_habitat" %in% colnames(species)){
                        species <- species |>
                          dplyr::rename(Source_for_predominant_land_use = Source_for_predominant_habitat)
                      }
                      
                      if ("Eco_region_distance_metres" %in% colnames(species)){
                        species <- species |>
                          dplyr::rename(Ecoregion_distance_metres = Eco_region_distance_metres)
                      }
                      if ("Max_linear_extent" %in% colnames(species)){
                        species <- species |>
                          dplyr::rename(Max_linear_extent_metres = Max_linear_extent)
                      }
                      
                      
                      # If available, load sites data
                      if (paths$sites != ""){
                        # ensure paths$sites is RDS
                        assertthat::assert_that(
                          assertthat::has_extension(paths$sites, "rds")
                          )
                      
                        sites <- readr::read_rds(paths$sites)
                        sites$version <- version
                        
                        if ("N_samples" %in% colnames(sites)){
                          if (inherits(sites$N_samples, "table")){
                            sites$N_samples <- as.vector(sites$N_samples)
                          }  
                        }
                        
                        assertthat::assert_that(
                          assertthat::has_name(sites, "Source_ID"),
                          assertthat::has_name(sites, "SS"),
                          assertthat::has_name(sites, "Diversity_metric"),
                          assertthat::has_name(sites, "SSBS"),
                          assertthat::has_name(sites, "Predominant_land_use") || assertthat::has_name(sites, "Predominant_habitat"),
                          assertthat::has_name(sites, "Use_intensity")
                        )
                        
                        # Ensure that Predominant_land_use is present - for LUI
                        if (!"Predominant_land_use" %in% colnames(sites)){
                          species <- species |>
                            dplyr::mutate(Predominant_land_use = Predominant_habitat)
                        }
                      }
                      
                    
                      # if site path missing, generate sites 
                      if (paths$sites == ""){
                        
                        assertthat::assert_that(
                          assertthat::has_name(species, "Sample_start_earliest"),
                          assertthat::has_name(species, "Sample_end_latest"),
                          assertthat::has_name(species, "Higher_taxon"),
                          assertthat::has_name(species, "Predominant_land_use"),
                          assertthat::has_name(species, "Use_intensity")
                          )
                        message(crayon::green("Generating sites dataframe...\n"))
                        
                        
                        # if sample_midpoint not present and sample_strt/end are...
                        if ((!"Sample_midpoint" %in% colnames(species) &
                             "Sample_start_earliest" %in% colnames(species) & 
                             "Sample_end_latest" %in% colnames(species))){
                          species <- species |>
                            dplyr::mutate(Sample_midpoint = Sample_start_earliest + floor((Sample_end_latest - Sample_start_earliest)/2))
                        }
                        
                        site_specific_cols <- c("Source_ID", "Reference",                          
                                                "Study_number", "Study_name",                             
                                                "SS", "Diversity_metric",                       
                                                "Diversity_metric_unit", 
                                                "Diversity_metric_type",                  
                                                "Diversity_metric_is_effort_sensitive",    
                                                "Diversity_metric_is_suitable_for_Chao",
                                                "Sampling_method","Sampling_effort_unit",                   
                                                "Study_common_taxon", 
                                                "Rank_of_study_common_taxon",
                                                "Site_number", "Site_name",                              
                                                "Block", "SSS",                                    
                                                "SSB", "SSBS",                                   
                                                "Sample_start_earliest", "Sample_end_latest",
                                                "Sample_midpoint",
                                                "Sample_date_resolution",                 
                                                "Max_linear_extent", "Habitat_patch_area_square_metres",
                                                "Sampling_effort", 
                                                # "Rescaled_sampling_effort",       
                                                "Habitat_as_described", 
                                                "Predominant_land_use",                   
                                                "Source_for_predominant_land_use",
                                                "Use_intensity",  
                                                # "Fragmentation_layout", 
                                                "Km_to_nearest_edge_of_habitat",
                                                "Years_since_fragmentation_or_conversion",
                                                "Transect_details", "Coordinates_method",                     
                                                "Longitude", "Latitude",                               
                                                "Country_distance_metres", "Country",                                
                                                "UN_subregion", "UN_region",                              
                                                "Ecoregion_distance_metres", "Ecoregion",                              
                                                "Biome", "Realm",                                  
                                                "Hotspot", "Wilderness_area")
                        
                        # subset site_specific cols to those in species??
                        site_specific_cols <- site_specific_cols[site_specific_cols %in% colnames(species)]
                        
                        # Summarise to site info
                        sites <- species |>
                          dplyr::group_by(dplyr::across(dplyr::all_of(site_specific_cols))) |>
                          dplyr::summarise(Higher_taxa = paste0(unique(Higher_taxon), collapse = ","),
                                           Sample_midpoint = Sample_midpoint[1], 
                                           N_samples = dplyr::n()) |>
                          dplyr::ungroup()
                        
    
                        sites$version <- version
                      }
                      
                      # %%%%%
                      # Add LUI column - Only if no LUI already
                      if (!"LUI" %in% colnames(sites)){
                        sites <- sites |>
                          dplyr::mutate(LUI = paste(Predominant_land_use, Use_intensity,
                                                    sep = "_"))
                      }
                      if (!"LUI" %in% colnames(species)){
                        species <- species |>
                          dplyr::mutate(LUI = paste(Predominant_land_use, Use_intensity,
                                                    sep = "_"))
                      }
                      
                      # %%%%%
                      
                      
                      # And check 
                      # Outsourced security checks
                      self$checkself(sites, species)
                      
                      # Save inputs in list
                      self$version <- version
                      self$last_update <- lubridate::today() # lubridate::now()
                      self$paths <- paths
                      self$sites <- sites
                      self$species <- species
                      
                      # Add some data processing flags
                      self$processing_flags <- list(measures_effort_corrected = F,
                                                    abundance_rescaled = F,
                                                    filtering_steps = list(),
                                                    joined_information = list(),
                                                    lu_legend_source = "Original LU/I legend."
                                                    )
                    },
                    
                    #' @description
                    #' Self check function to ensure that the site and species
                    #' data.frames contain the correct columns and align in content.
                    #'
                    #' @param sites (`data.frame`).
                    #' @param species (`data.frame`).
                    #' @return (`NULL`).
                    checkself = function(sites = self$sites, species = self$species){
                      # ------ #
                      # Security checks
                      assertthat::assert_that(
                        assertthat::has_name(sites, "Source_ID"),
                        assertthat::has_name(sites, "SS"),
                        assertthat::has_name(sites, "Diversity_metric"),
                        assertthat::has_name(sites, "SSBS"),
                        assertthat::has_name(sites, "Predominant_land_use") || assertthat::has_name(sites, "Predominant_habitat"),
                        assertthat::has_name(sites, "Use_intensity"),
                        # TODO: extra essential columns
                        assertthat::has_name(sites, "SSB"), 
                        assertthat::has_name(sites, "LUI"), 
                        assertthat::has_name(sites, "Longitude"), 
                        assertthat::has_name(sites, "Latitude"),

                        # And the species data
                        assertthat::has_name(species, "Source_ID"),
                        assertthat::has_name(species, "SS"),
                        assertthat::has_name(species, "SSBS"),
                        assertthat::has_name(species, "Best_guess_binomial"),
                        assertthat::has_name(species, "Measurement"),
                        # Predominant_land_use, Use_intensity, LUI - optional, if present in sites,
                        # can be added from there... - could even be added later... via set_lu_legend
                        
                        assertthat::has_name(species, "Predominant_land_use"), 
                        assertthat::has_name(species, "Use_intensity"), 
                        assertthat::has_name(species, "LUI"), 

                        assertthat::has_name(species, "Diversity_metric_type"), 

                        # Should cols needed for div metrics be enforced here...
                        # What use is the database if metrics can't be calculated..?!
                        # Taxon_name_entered, Best_guess_binomial, Sampling_effort
                        assertthat::has_name(species, "Taxon_name_entered"), 
                        # assertthat::has_name(species, "Best_guess_binomial"), 
                        assertthat::has_name(species, "Sampling_effort")
                        # Diversity_metric_is_effort_sensitive

                      )
                      # Finally, check that all sites are in species
                      # NOTE RC - source ID is just the paper...
                      assertthat::assert_that(
                        nrow(sites) > 0, nrow(species) > 0,
                        all( unique(sites$Source_ID) %in% unique(species$Source_ID) )
                      )
                      invisible()
                    },
                    
                    #' @description
                    #' Base printout function presents core information about the 
                    #' Database object.
                    #'
                    #' @return (`NULL`).
                    print = function() {
                      # Compute some stats
                      
                      # unique(colnames(self$site_contrasts)) 
                      # Get rid of basic details, check for columns related to metrics...
                      
                      site_contrasts_poss <- c("LUI_contrast", 
                                               "Geogr_dist_hav", "Geogr_dist_geod",
                                               "Latitude_midpoint", "Longitude_midpoint")
                      
                      stats <- list(number_sites = nrow(self$sites),
                                    number_lu = length( unique(self$sites$Predominant_land_use) ),
                                    number_lui = length(unique(self$sites$Use_intensity)),
                                    number_species = length( unique(self$species$Best_guess_binomial) ),
                                    filters_applied = paste(self$processing_flags$filtering_steps, collapse = "; "),
                                    # div_metrics = paste(sort(unique(self$diversity$Metric)), collapse = ", ") 
                                    
                                    joined_information = paste(self$processing_flags$joined_information, collapse = "; "),
                                    lu_legend_source = self$processing_flags$lu_legend_source,

                                    div_within_metrics = paste(sort(unique(self$diversity_within$Metric)), collapse = ", "),
                                    div_between_metrics = paste(sort(unique(self$diversity_between$Metric)), collapse = ", "),
                                    # TODO: site contrast info
                                    site_contrast_metrics = paste(sort(unique(colnames(self$site_contrasts)[colnames(self$site_contrasts) %in% site_contrasts_poss])), 
                                                                  collapse = ", ")
                                    )  
                      
                      
                      # Neaten diversity metrics formatting
                      if (is.null(self$diversity_within)){
                        stats$div_within_metrics <- "NA"
                      }
                      if (!is.null(self$diversity_within)){
                        stats$div_within_metrics <- gsub("_", " ", stats$div_within_metrics)
                        stats$div_within_metrics <- gsub("Simpson diversity", "Simpson\'s diversity",
                                                         stats$div_within_metrics)
                      }
                      
                      if (is.null(self$diversity_between)){
                        stats$div_between_metrics <- "NA"
                      }
                      if (!is.null(self$diversity_between)){
                        stats$div_between_metrics <- gsub("_", " ", stats$div_between_metrics)
                        # stats$div_between_metrics <- gsub("Simpson diversity", "Simpson\'s diversity",
                        #                           stats$div_between_metrics)
                      }
                      
                      if (is.null(self$site_contrasts)){
                        stats$site_contrast_metrics <- "NA"
                      }
                      if (!is.null(self$site_contrasts)){
                        # stats$site_contrast_metrics <- gsub("lui_contrast", "LUI contrast",
                        #                                  stats$site_contrast_metrics)
                        stats$site_contrast_metrics <- gsub("Geogr_dist_hav", "Geographic dist. (haversine)",
                                                            stats$site_contrast_metrics)
                        stats$site_contrast_metrics <- gsub("Geogr_dist_geod", "Geographic dist. (geodesic)",
                                                            stats$site_contrast_metrics)
                        stats$site_contrast_metrics <- gsub("_", " ",
                                                            stats$site_contrast_metrics)
                      }
                      
                      # Neaten filtering formatting
                      if (length(self$processing_flags$filtering_steps) == 0){
                        stats$filters_applied <- "NA"
                      }

                      if (length(self$processing_flags$joined_information) == 0){
                        stats$joined_information <- "NA"
                      }
                      message(
                        paste0("Hibiiscus PREDICTS-like database (including version/s: ", crayon::bold(self$version),")","\n",
                               "(Last updated: ", crayon::green( self$last_update ), ")",
                               "\n",
                               "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                               "\t Number of species:             ", crayon::green(stats$number_species), "\n",

                               "\t Land use (class x intensity):  ", crayon::green(stats$number_lu), " x ",
                                                                     crayon::green(stats$number_lui), "\n",

                               "\t LU/I legend sourced from:      ", crayon::green(stats$lu_legend_source), "\n",
                               
                               "\n",
                               
                               "\t Effort-correction applied:     ", crayon::green(self$processing_flags$measures_effort_corrected), "\n",
                               
                               "\t Filtering conducted:           ", crayon::green(gsub("; ", 
                                                                                        ";\n\t                                ",
                                                                                        stats$filters_applied)), "\n",

                               "\t Data added to database:        ", crayon::green(gsub("; ", 
                                                                                        ";\n\t                                ",
                                                                                        stats$joined_information)), "\n",

                               "\n",

                               "\t Within-site diversity metrics calculated:   ", crayon::green(stats$div_within_metrics), "\n",
                               "\t Between-site diversity metrics calculated:  ", crayon::green(stats$div_between_metrics), "\n",
                               "\t Site-site contrast metrics calculated:      ", crayon::green(stats$site_contrast_metrics), "\n"
                               )
                        )
                      
                    },
                    
                    #' @description
                    #' Function to save the Database object to the specified 
                    #' file-path ('output_path', must be '.rds').
                    #'
                    #' @param output_path (`character()`)
                    #' @return (`NULL`).
                    save_hib = function(output_path = NULL){
                      # hib_save ??
                      assertthat::assert_that(
                        is.character(output_path)
                      )
                      # if (is.null(output_path)){
                      #   output_path <- self$paths$path_output
                      # }
                      
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
                      
                      message(crayon::green(paste0("Saving hib database to: '", output_path, "'\n")) )
                      
                      # Note: currently doesn't make the folder or check parent folders exist...
                      # Could include - with a warning...

                      saveRDS(self, output_path)
                      assertthat::assert_that(file.exists(output_path),
                                              msg = "Writing a copy of the database failed...")
                      invisible()
                    },
                    
                    #' @description
                    #' Function to apply effort correction to the species data. 
                    #' Uses a simple, linear scaling, such that sampling effort
                    #' of a site is scaled relative to the maximum value per study,
                    #' and multiplied by the original measure of diversity.
                    #'
                    #' @return (`HibDatabase`).
                    hib_effort_correction = function(){
                      # Assume data is in self$species?
                      # Check that relevant columns are present
                      assertthat::assert_that(
                        assertthat::has_name(self$species, "Sampling_effort"),
                        assertthat::has_name(self$species, "SS"),
                        assertthat::has_name(self$species, "Measurement"),
                        assertthat::has_name(self$species, "Diversity_metric_is_effort_sensitive")
                      )

                      # add a message...
                      message(crayon::green("Correcting measuremnts for sampling effort...\n"))
                      
                      self$species <- self$species |>
                        dplyr::mutate(Sampling_effort = dplyr::case_when(is.na(Sampling_effort) ~ 1,
                                                                         T ~ Sampling_effort)) |>
                        dplyr::group_by(SS) |>
                        dplyr::mutate(Rescaled_sampling_effort = Sampling_effort / max(Sampling_effort)) |>
                        dplyr::ungroup() |>
                        dplyr::mutate(Effort_corrected_measurement = dplyr::case_when(Diversity_metric_is_effort_sensitive == T ~ Measurement / Rescaled_sampling_effort,
                                                                                      T ~ Measurement))
                      # Add a flag to some background info - Measures are effort corrected = T
                      self$processing_flags$measures_effort_corrected <- T
                      self$last_update <- lubridate::today()
                      invisible(self)
                    } ,
                    
                    #' @description
                    #' Function to generate figures summarising the data present
                    #' in the HibDatabase.
                    #' 
                    #' @param type (`character()`) Options are: "map-sites", 
                    #' "map-studies", "biome-N", "ecoregion-N", "realm-N", "class-N".
                    #' @return (`ggplot`).
                    plot = function(type = "map-sites"){
                      # Figure options
                      # - site map (map-sites)
                      # - study map (map-studies)
                      # - country heatmap - N sites (heatmap-sites-countries)
                      # - country heatmap - N studies (heatmap-studies-countries)
                      # - biome N (biome-N)
                      # - ecoregion N (ecoregion-N)
                      # - continent N - map/bar? (continent-N) ??
                      # - realm N (realm-N)
                      # - class N (class-N)
                      
                      # TODO: what about within-site and between-site diversity metrics??
                      
                      # Assert that type arg matches one of...
                      type <- match.arg(type,
                                        c("map-sites", "map-studies", 
                                          "biome-N", "ecoregion-N", "realm-N",
                                          "class-N"),
                                        several.ok = F)
                      # Clr blind friendly...
                      # black, orange, sky blue, bluish green, yellow, blue, vermillion, reddish purple
                      # pal <- c("#000000", "#e69f00", "#56b4e9", "#009e73",
                      #                   "#f0e442", "#0072b2", "#d55e00", "#cc79a7")
                                        
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
                      
                      # TODO: combine similar plot types into one call
                      # %%%%%
                      if (type %in% c("map-sites", "map-studies")){

                        # Assert that Lat and long are in sites
                        assertthat::assert_that(
                          assertthat::has_name(self$sites, "Latitude"),
                          assertthat::has_name(self$sites, "Longitude")
                          )
                        # Just get a coarse-ish world map...
                        tmp <- rnaturalearth::ne_countries(scale = 50, returnclass = 'sf')
                        
                        if (type == "map-sites"){
                          plt_df <- self$sites
                        }
                        if (type == "map-studies"){
                          plt_df <- self$sites |>
                            dplyr::group_by(SS, version) |>
                            dplyr::summarise(N_sites = dplyr::n(),
                                             Longitude = mean(Longitude, na.rm = T),
                                             Latitude = mean(Latitude, na.rm = T))
                        }
                        # Get x/y limits for sites...
                        x_rng <- range(plt_df$Longitude, na.rm = T)
                        y_rng <- range(plt_df$Latitude, na.rm = T)
                        
                        # buffer by 5 degrees?
                        x_rng[1] <- round(x_rng[1] - 5, 2)
                        x_rng[2] <- round(x_rng[2] + 5, 2)
                        
                        y_rng[1] <- round(y_rng[1] - 5, 2)
                        y_rng[2] <- round(y_rng[2] + 5, 2)
                        
                        x_rng[1] <- max(c(-180, x_rng[1]))
                        x_rng[2] <- min(c(180, x_rng[2]))
                        
                        y_rng[1] <- max(c(-90, y_rng[1]))
                        y_rng[2] <- min(c(90, y_rng[2]))
                        
                        plt <- ggplot2::ggplot() +
                          ggplot2::geom_sf(data = tmp,
                                  fill = "lightgrey",
                                  alpha = 0.4,
                                  show.legend = F) +
                          ggplot2::xlim(x_rng) +
                          ggplot2::ylim(y_rng) +
                          ggplot2::labs(colour = "Dataset", #"PREDICTS\nversion",
                               x = ggplot2::element_blank(),
                               y = ggplot2::element_blank()) +
                          ggplot2::guides(colour = ggplot2::guide_legend(
                            override.aes=list(size = 3,
                                              alpha = 0.9))) +
                          ggplot2::theme_bw() +
                          basic_thm 
                        
                        if (type == "map-sites"){
                          plt <- plt +
                            ggplot2::geom_point(ggplot2::aes(x = Longitude, y  = Latitude, 
                                           colour = version),
                                       size = 0.15,
                                       alpha = 0.25,
                                       data = plt_df)
                          
                        }
                        if (type == "map-studies"){
                          plt <- plt +
                            ggplot2::geom_point(ggplot2::aes(x = Longitude, y  = Latitude,
                                           size = N_sites,
                                           colour = version),
                                       alpha = 0.25,
                                       data = plt_df) +
                            ggplot2::labs(size = "No. studies")
                        }
                      }
                      # %%%%%
                                            
                      # Issue with country-level heatmap - names in predicts data
                      # might not match those in natural earth data (admin/sovereignt)
                      # - country heatmap - N sites (heatmap-sites-countries)
                      # - country heatmap - N studies (heatmap-studies-countries)
                      
                      # %%%%%
                      if (type %in% c("biome-N", "ecoregion-N", "realm-N")){
                        var_ <- tools::toTitleCase(type) |> gsub("-N", "", .)

                        # Assert that Lat and long are in sites
                        assertthat::assert_that(
                          assertthat::has_name(self$sites, var_)
                          )
                        
                        count_df <- self$sites |>
                          dplyr::group_by(!!sym(var_), version) |>
                          dplyr::summarise(Sites = dplyr::n_distinct(SS),
                                           Studies = dplyr::n_distinct(SSBS)) |>
                          tidyr::pivot_longer(c(Sites, Studies),
                                              names_to = "dat_type",
                                              values_to = "N")
                        
                        plt <- ggplot2::ggplot(count_df) +
                          ggplot2::geom_bar(ggplot2::aes(x = !!sym(var_), y = N,
                                       fill = version),
                                   # fill = "#009E73",
                                   stat = "identity") +
                          ggplot2::labs(x = var_,
                               y = "N",
                               fill = "Dataset" #"PREDICTS\nversion"
                               ) +
                          ggplot2::facet_wrap(~dat_type,
                                     scales = "free_y") +
                          ggplot2::theme_bw() +
                          basic_thm +
                          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, 
                                                           hjust = 1, 
                                                           vjust = 1))
                        
                        if (length(unique(count_df[[var_]])) > 20){
                          plt <- plt +
                            ggplot2::theme(axis.text.x = ggplot2::element_blank())
                        }
                      }
                      # %%%%%
                      
                      # - class N (class-N)
                      if (type == "class-N"){
                        # Assert that Lat and long are in sites
                        assertthat::assert_that(
                          assertthat::has_name(self$sites, "Class")
                          )
                        class_df <- self$species |>
                          dplyr::group_by(Class, version) |>
                          dplyr::summarise(N_ = dplyr::n()) |>
                          # dplyr::filter(N_ > 0) |>
                          dplyr::filter(Class != "")
                        
                        plt <- ggplot2::ggplot(class_df) +
                          ggplot2::geom_bar(ggplot2::aes(x = Class, y = N_,
                                       fill = version),
                                   # fill = "#E69F00",
                                   stat = "identity") +
                          ggplot2::labs(x = "Class",
                               y = "No. records",
                               fill = "Dataset" #"PREDICTS\nversion"
                               ) +
                          ggplot2::theme_bw() +
                          basic_thm +
                          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, 
                                                           hjust = 1, 
                                                           vjust = 1))
                        
                        if (length(unique(class_df$Class)) > 20){
                          plt <- plt +
                            ggplot2::theme(axis.text.x = ggplot2::element_blank())
                        }
                      }
                      (plt)
                    }
                    # Filter database
                    # TODO: put filter methods here...
                    # by sites: "Country", "UN_subregion" , "UN_region",  "Ecoregion" , "Biome", "Realm"...
                    # hib_filter_sites()
                    
                  ),
                  # Private entries to be created
                  private = list(
                    siteids = list(),
                    length = function() base::length(self$sites$SS)
                  )
)

