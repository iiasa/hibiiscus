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
                    # diversity = NULL,
                    diversity_within = NULL,
                    diversity_between = NULL,
                    
                    site_contrasts = NULL,
                    
                    processing_flags = NULL,
                    # Executed when created
                    initialize = function(version, paths){
                      assertthat::assert_that(is.character(version), is.list(paths))
                      
                      # Assert path$sites/species is character ?? - already checked in create_hib
    
                      # ------ #
                      message(crayon::green(paste0("Initializing input datasets for PREDICTS version ", version, "...")) )
                      
                      # Now load both datasets
                      # Note RC: But the species data has all the sites
                      # And the 2022 data is just the full spp data
                      
            
                      # %%%%%
                      # Load species data
                      species <- readr::read_rds(paths$species)
                      species$version <- version
                      
                      # If predicts2022 column names, rename to align with 2016...
                      if ("Predominant_habitat" %in% colnames(species)){
                        species <- species %>%
                          dplyr::rename(Predominant_land_use = Predominant_habitat)
                      }
                      if ("Source_for_predominant_habitat" %in% colnames(species)){
                        species <- species %>%
                          dplyr::rename(Source_for_predominant_land_use = Source_for_predominant_habitat)
                      }
                      
                      if ("Eco_region_distance_metres" %in% colnames(species)){
                        species <- species %>%
                          dplyr::rename(Ecoregion_distance_metres = Eco_region_distance_metres)
                      }
                      if ("Max_linear_extent" %in% colnames(species)){
                        species <- species %>%
                          dplyr::rename(Max_linear_extent_metres = Max_linear_extent)
                      }
                      
                      
                      # If available, load sites data
                      if (paths$sites != ""){
                        sites <- readr::read_rds(paths$sites)
                        sites$version <- version
                        
                        if ("N_samples" %in% colnames(sites)){
                          if (class(sites$N_samples) == "table"){
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
                      }
                      
                    
                      # if site path missing, generate sites 
                      if (paths$sites == ""){
                        
                        assertthat::assert_that(
                          assertthat::has_name(species, "Sample_start_earliest"),
                          assertthat::has_name(species, "Sample_end_latest"),
                          assertthat::has_name(species, "Higher_taxon"),
                          assertthat::has_name(species, "Predominant_land_use"),
                          assertthat::has_name(species, "Use_intensity")
                          # Predominant land use
                          # Use intensity
                          )
                        message(crayon::green("Generating sites dataframe..."))
                        
                        
                        # if sample_midpoint not present and sample_strt/end are...
                        if ((!"Sample_midpoint" %in% colnames(species) &
                             "Sample_start_earliest" %in% colnames(species) & 
                             "Sample_end_latest" %in% colnames(species))){
                          species <- species %>%
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
                        sites <- species %>%
                          dplyr::group_by(dplyr::across(dplyr::all_of(site_specific_cols))) %>%
                          dplyr::summarise(Higher_taxa = paste0(unique(Higher_taxon), collapse = ","),
                                           Sample_midpoint = Sample_midpoint[1], 
                                           N_samples = n()) %>%
                          dplyr::ungroup()
                        
    
                        sites$version <- version
                      }
                      
                      
                      # Add LUI column - 
                      sites <- sites %>%
                        mutate(LUI = paste(Predominant_land_use, Use_intensity,
                                           sep = "_"))
                      species <- species %>%
                        mutate(LUI = paste(Predominant_land_use, Use_intensity,
                                           sep = "_"))
                      
                      # Alternative initialisation approach - version specific
                      # # Use some standard site columns - subset to those in species df
                      # 
                      # # Check for data structure... - if core columns missing break
                      # # - e.g. Source_ID, SS, SSBS, Measure, Predominant_land_use, 
                      # #         Diversity_metric, Diversity_metric_type, Best_guess_binomial
                      # 
                      # 
                      # # Add version value to dfs
                      # 
                      # # %%%%%
                      # 
                      # 
                      # if(version == "2016"){
                      #   sites <- readr::read_rds(paths$sites)
                      #   species <- readr::read_rds(paths$species)
                      #   # Add version info to sites/species
                      # }
                      # if(version == "2022"){
                      #   
                      #   species <- readr::read_rds(paths$species)
                      #   
                      #   message(crayon::green("Generating sites dataframe..."))
                      #   
                      #   # unique values of important columns...
                      #   # also need to generate N_samples, Sample_midpoint, Higher_taxa
                      #   # and maybe Rescaled_sampling_effort
                      #   
                      #   # TODO: neater way to do this?? - all cols except 
                      #   # -Higher_taxon
                      #   # -...
                      #   site_specific_cols <- c("Source_ID", "Reference",                          
                      #                           "Study_number", "Study_name",                             
                      #                           "SS", "Diversity_metric",                       
                      #                           "Diversity_metric_unit", 
                      #                           "Diversity_metric_type",                  
                      #                           "Diversity_metric_is_effort_sensitive",    
                      #                           "Diversity_metric_is_suitable_for_Chao",
                      #                           "Sampling_method","Sampling_effort_unit",                   
                      #                           "Study_common_taxon", 
                      #                           "Rank_of_study_common_taxon",
                      #                           "Site_number", "Site_name",                              
                      #                           "Block", "SSS",                                    
                      #                           "SSB", "SSBS",                                   
                      #                           "Sample_start_earliest", "Sample_end_latest",
                      #                           "Sample_date_resolution",                 
                      #                           "Max_linear_extent", "Habitat_patch_area_square_metres",
                      #                           "Sampling_effort", 
                      #                           # "Rescaled_sampling_effort",       
                      #                           "Habitat_as_described", 
                      #                           "Predominant_habitat",                   
                      #                           "Source_for_predominant_habitat",
                      #                           "Use_intensity",                          
                      #                           "Km_to_nearest_edge_of_habitat",
                      #                           "Years_since_fragmentation_or_conversion",
                      #                           "Transect_details", "Coordinates_method",                     
                      #                           "Longitude", "Latitude",                               
                      #                           "Country_distance_metres", "Country",                                
                      #                           "UN_subregion", "UN_region",                              
                      #                           "Eco_region_distance_metres", "Ecoregion",                              
                      #                           "Biome", "Realm",                                  
                      #                           "Hotspot", "Wilderness_area")
                      #   
                      #   # subset site_specific cols to those in species??
                      #   
                      #   # if sample_midpoint not present and sample_strt/end are...
                      #   sites <- species %>%
                      #     dplyr::mutate(Sample_midpoint = Sample_start_earliest + floor((Sample_end_latest - Sample_start_earliest)/2)) %>%
                      #     dplyr::group_by(dplyr::across(dplyr::all_of(site_specific_cols))) %>%
                      #     dplyr::summarise(Higher_taxa = paste0(unique(Higher_taxon), collapse = ","),
                      #                      Sample_midpoint = Sample_midpoint[1], 
                      #                      N_samples = n()) %>%
                      #     dplyr::ungroup()
                      #   
                      #   # Name changes to match PREDICTS 2016
                      #   # Convert Predominant_habitat to Predominant_land_use
                      #   # Convert Eco_region_distance_metres to Ecoregion_distance_metres
                      #   # Convert Max_linear_extent to Max_linear_extent_metres
                      #   species <- species %>%
                      #     dplyr::rename(Predominant_land_use = Predominant_habitat,
                      #                   Ecoregion_distance_metres = Eco_region_distance_metres,
                      #                   Max_linear_extent_metres = Max_linear_extent)
                      #   sites <- sites %>%
                      #     dplyr::rename(Predominant_land_use = Predominant_habitat,
                      #                   Ecoregion_distance_metres = Eco_region_distance_metres,
                      #                   Max_linear_extent_metres = Max_linear_extent)
                      #   
                      #   # Differences in column names
                      #   # 2022 version v 2016 version
                      #   # Predominant_habitat (2022) is equiv to Predominant_land_use 
                      #   # - but Predominant_habitat splits primary forest and primary non-forest
                      #   # 2022 has 11 levels, 2016 has 10
                      #   # Can split primary into forest/non-forest based on biome ??
                      #   
                      #   # Source_for_predominant_habitat =? Source_for_predominant_land_use
                      #   
                      #   # Max_linear_extent is in 2022, Max_linear_extent_metres is in 2016...
                      #   
                      #   # Eco_region_distance_metres =? Ecoregion_distance_metres
                      #   
                      #   # Higher_taxa -- "," separated character of higher taxa present in site
                      #   
                      #   # Sample_midpoint missing from 2022 - sample_end + sample_strt / 2 ??
                      #   
                      #   # Rescaled_sampling_effort and N_sampled missing 
                      #   # Rescaled_sampling_effort can be easily calculated - needed?
                      #   
                      #   
                      #   # Add version info to sites/species
                      # }
                      # 
                      # # TODO:
                      # # Add custom data source version - but retain need for sites and species?
                      # # Definitely need species.., sites is optional...
                      # if(version == "custom"){
                      #   sites <- readr::read_rds(paths$sites)
                      #   species <- readr::read_rds(paths$species)
                      #   # Add version info to sites/species
                      # }
                      
                      # And check 
                      # Outsourced security checks
                      self$checkself(sites, species)
                      
                      # Save inputs in list
                      self$version <- version
                      self$last_update <- lubridate::today()
                      self$paths <- paths
                      self$sites <- sites
                      self$species <- species
                      
                      # Add some data processing flags
                      self$processing_flags <- list(measures_effort_corrected = F,
                                                    abundance_rescaled = F,
                                                    filtering_steps = list())
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
                        assertthat::has_name(sites, "Predominant_land_use") || assertthat::has_name(sites, "Predominant_habitat"),
                        assertthat::has_name(sites, "Use_intensity"),
                        
                        # And the species data
                        assertthat::has_name(species, "Source_ID"),
                        assertthat::has_name(species, "SS"),
                        assertthat::has_name(species, "SSBS"),
                        assertthat::has_name(species, "Best_guess_binomial"),
                        assertthat::has_name(species, "Measurement")
                      )
                      # Finally, check that all sites are in species
                      # NOTE RC - source ID is just the paper...
                      assertthat::assert_that(
                        nrow(sites)>0, nrow(species)>0,
                        all( unique(sites$Source_ID) %in% unique(species$Source_ID) )
                      )
                      invisible()
                    },
                    # Base printout function
                    print = function() {
                      # Compute some stats
                      
                      # unique(colnames(self$site_contrasts)) 
                      # Get rid of basic details, check for columns related to metrics...
                      
                      stats <- list(number_sites = nrow(self$sites),
                                    number_lu = length( unique(self$sites$Predominant_land_use) ),
                                    number_lui = length(unique(self$sites$Use_intensity)),
                                    number_species = length( unique(self$species$Best_guess_binomial) ),
                                    filters_applied = paste(self$processing_flags$filtering_steps, collapse = "; "),
                                    # div_metrics = paste(sort(unique(self$diversity$Metric)), collapse = ", ") 
                                    
                                    div_within_metrics = paste(sort(unique(self$diversity_within$Metric)), collapse = ", "),
                                    div_between_metrics = paste(sort(unique(self$diversity_between$Metric)), collapse = ", ")
                                    # TODO: site contrast info
                                    # site_contrast_metrics = unique(colnames(self$site_contrasts))
                                    )  
                      # TODO: check both within and between site metrics/diversity data
                      
                      # Neaten diversity metrics formatting
                      if (is.null(self$diversity_within)){
                        stats$div_within_metrics = "NA"
                      }
                      if (!is.null(self$diversity_within)){
                        stats$div_within_metrics <- gsub("_", " ", stats$div_within_metrics)
                        stats$div_within_metrics <- gsub("Simpson diversity", "Simpson\'s diversity",
                                                         stats$div_within_metrics)
                      }
                      
                      if (is.null(self$diversity_between)){
                        stats$div_between_metrics = "NA"
                      }
                      if (!is.null(self$diversity_between)){
                        stats$div_between_metrics <- gsub("_", " ", stats$div_between_metrics)
                        # stats$div_between_metrics <- gsub("Simpson diversity", "Simpson\'s diversity",
                        #                           stats$div_between_metrics)
                      }
                      
                      
                      # Neaten filtering formatting
                      if (length(self$processing_flags$filtering_steps) == 0){
                        stats$filters_applied = "NA"
                      }
                      message(
                        paste0("Hibiiscus PREDICTS database (including versions: ", crayon::bold(self$version),")","\n",
                               "(Last updated: ", crayon::green( self$last_update ), ")",
                               "\n",
                               "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                               "\t Land use (class x intensity):  ", crayon::green(stats$number_lu), " x ",
                                                                     crayon::green(stats$number_lui), "\n",
                               "\t Number of species:             ", crayon::green(stats$number_species), "\n",
                               "\n",
                               
                               "\t Effort-correction applied:     ", crayon::green(self$processing_flags$measures_effort_corrected), "\n",
                               
                               "\t Filtering conducted:           ", crayon::green(gsub("; ", 
                                                                                        ";\n\t                                ",
                                                                                        stats$filters_applied)), "\n",
                               "\t Within-site diversity metrics calculated:  ", crayon::green(stats$div_within_metrics), "\n",
                               "\t Between-site diversity metrics calculated:  ", crayon::green(stats$div_between_metrics), "\n"
                               )
                        )
                      
                    },
                    # Save object
                    save_hib = function(output_path = NULL){
                      assertthat::assert_that(
                        is.character(output_path)
                      )
                      # if (is.null(output_path)){
                      #   output_path <- self$paths$path_output
                      # }
                      
                      if(!assertthat::has_extension(output_path, "rds")){
                        output_path <- paste0(output_path, ".rds")
                      }
                      
                      if (file.exists(output_path)){
                        message(crayon::yellow(paste0("Warning: '", output_path, "' already exists.",
                                                      "\n",
                                                      "This file will be updated.\n")))
                      }
                      
                      message(crayon::green(paste0("Saving hib database to: '", output_path, "' ...\n")) )
                      
                      saveRDS(self, output_path)
                      assertthat::assert_that(file.exists(output_path),
                                              msg = "Writing a copy of the database failed...")
                      invisible()
                    },
                    # Correct measurement for sampling effort
                    hib_effort_correction = function(){
                      # Assume data is in self$species?
                      # Check that relevant columns are present
                      assertthat::assert_that(
                        assertthat::has_name(self$species, "Sampling_effort"),
                        assertthat::has_name(self$species, "SS"),
                        assertthat::has_name(self$species, "Measurement")
                      )

                      # add a message...
                      message(crayon::green("Correcting measuremnts for sampling effort...\n"))
                      
                      self$species <- self$species %>%
                        dplyr::mutate(Sampling_effort = dplyr::case_when(is.na(Sampling_effort) ~ 1,
                                                                         T ~ Sampling_effort)) %>%
                        dplyr::group_by(SS) %>%
                        dplyr::mutate(Rescaled_sampling_effort = Sampling_effort / max(Sampling_effort)) %>%
                        dplyr::ungroup() %>%
                        dplyr::mutate(Effort_corrected_measurement = dplyr::case_when(Diversity_metric_is_effort_sensitive == T ~ Measurement / Rescaled_sampling_effort,
                                                                                      T ~ Measurement))
                      # Add a flag to some background info - Measures are effort corrected = T
                      self$processing_flags$measures_effort_corrected <- T
                      self$last_update <- lubridate::today()
                      invisible(self)
                    } #,
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

