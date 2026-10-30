#' create_dummy_data: Internal function to create dummy data for PREDICTS-style models
#' 
#' @description
#' The package itself should not contain or reference any external datasets
#' or their paths (which are system specific anyway).
#' Rather a dummy dataset with the same format should be used, also allowing 
#' us to create tests for the package that are not reliant on external data.
#' 
#' @details
#' This function creates dummy data in the format of the PREDICTS database. 
#' A separate function (TBD) can be used to create external information for 
#' projection or integration
#' 
#' @param what [`character`] stating what exactly to create. Available options are
#' \code{"sites"} to emulate a PREDICTS site database and \code{"species"} for a 
#' species database.
#' @param n_studies [`numeric`] on the number of studies to be created.
#' @param n_sites [`numeric`] Number of sites per study to be created.
#' @param lu_legend [`character`] of the land use legend to be generated for the 
#' dummy data.
#' @param ui_legend [`character`] of the use-intensity legend to be generated for the 
#' dummy data.
#' 
#' @return data.frame
#' 
#' @examples
#' \dontrun{
#'   out <- create_dummy_data()
#'   # Save in hibiscus data folder (lazy loading enabled)
#'   save(out, "data/dummy_data.rda")
#' }
#' 
#' @keywords internal
#' 
#' @author Martin Jung
#' @noRd
create_dummy_data <- function(what = "sites",
                              n_studies = 20, 
                              n_sites   = 10,
                              lu_legend = c("Primary", "Secondary", "Cropland", "Urban"),
                              ui_legend = c("Minimal_use", "Light_use", "Intense_use")
                              # TODO:
                              # lu_effects = c(1, 0.8, 0.5, 0.6),
                              # ui_effects = c(1, 0.75, 0.4)
                              # UN_region?
                              # study_radius
                              ){
  # Check input parameters
  assertthat::assert_that(
    is.character(what),
    is.numeric(n_studies), n_studies >= 3,
    is.numeric(n_sites), n_sites >= 2,
    is.character(lu_legend), length(lu_legend) >= 2,
    is.character(ui_legend), length(ui_legend) >= 1
  )
  # Match to available options
  what <- match.arg(what, c("sites", "species"), several.ok = FALSE)
  
  # print(what)
  # Make the land-use legend a factor, with the first level being the reference
  lu_legend <- factor(lu_legend, levels = c(lu_legend[1], lu_legend[2:length(lu_legend)]))
  
  # Normalize function
  norm <- function(y) (y - min(y)) / (max(y) - min(y)) 
  
  # Create a new data.frame container
  if(what == "sites"){
    
    # Create ids of studies
    ns <- vapply(seq_len(n_studies), function(study_index) {
      paste0(sample(c(letters, LETTERS, 0:9), 7L, replace = TRUE), collapse = "")
    }, character(1))
    assertthat::assert_that(length(unique(ns)) == n_studies)
    # lui <- c("Minimal use", "Light use", "Intense use")
   
    out <- tidyr::crossing(
      SS = paste0(ns, " ", 1),
      SSBS = 1:n_sites
    ) |> 
      dplyr::mutate(SSBS = paste0(SS, " ", SSBS),
                    # TODO: Add blocking scheme
                    Source_ID = substr(SS, start = 0, stop = 7),
                    Predominant_land_use = sample(lu_legend, size = n_sites * n_studies, replace = TRUE),
                    Use_intensity = sample(ui_legend, size = n_sites * n_studies, replace = TRUE, prob = c(1,.5,.2))
                    )

    # Other study settings
    # Normalize and invert
    out$Rescaled_sampling_effort <- 1- norm(stats::rexp(nrow(out),rate = 2))
    # NOTE: but minimum rescaled sample effort will currently be 0... - issue for 
    # rescaling abundance metrics
    # Alternative:
    # out$Sampling_effort <- rexp(nrow(out),rate = 2)
    # Rescaled sampling effort calculated within processing - not needed here
    # out <- out |>
    #   dplyr::group_by(SS) |>
    #   dplyr::mutate(Rescaled_sampling_effort = Sampling_effort / max(Sampling_effort)) |>
    #   dplyr::ungroup()
    
    
    # Generate dummy coordinates
    # TODO:
    # Random sample from terrestrial locations - per SS
    # Then randomly sample around that.. within some buffered distance
    
    ## Generate some random measurements ##
    out$Species_richness <- stats::rpois(n = nrow(out), 5) 
    out$LogAbundance <- stats::rnorm(n = nrow(out), 10, .5)
     
    # Just use abundance. - log-normal
    # out$Diversity_metric_type <- "Abundance"
    # out$Measurement <- rnorm(n = nrow(out),10, .5)
    
    # NOTE: should just make a species dataframe - then sites can be created by aggregating
    # generate columns for - 
    # Source_ID, SS, SSB, SSBS, Latitude, Longitude, 
    # Predominant_land_use, Use_intensity, LUI, 
    # Sampling_effort, Diversity_metric, Diversity_metric_type, Measurement,
    # Taxon_name_entered, Best_guess_binomial
    # Sample_start_earliest, Sample_end_latest, Higher_taxon
    
    
    # %%%%%%%%%%
    # TODO: hib Database class expects certain columns, and others are required downstream too!!
    
    # Sites: 
    # Source_ID, SS, SSBS, Diversity_metric, Predominant_land_use, Use_intensity
    
    # For div metrics: LUI
    # For modelling: SSB
    # For INLA: "Longitude", "Latitude"
    
    # Species: 
    # Source_ID, SS, Best_guess_binomial, Measurement
    
    # If sites is not provided, species also needs:
    # Sample_start_earliest, Sample_end_latest, Higher_taxon, Predominant_land_use, Use_intensity
    
    # For diversity metrics calc, species also needs:
    # SSBS, Diversity_metric_type - within_site
    # LUI, Taxon_name_entered - between_site
    # Measurement
    
    # For modelling: Best_guess_binomial, Sampling_effort, Taxon_name_entered
    
    # %%%%%%%%%%

  } else {
    # %%%%%%%%%%%%%%%%
    # for dev testing
    # n_studies <- 20
    # n_sites = 10
    # lu_legend = c("Primary", "Secondary", "Cropland", "Urban")
    # ui_legend = c("Minimal_use", "Light_use", "Intense_use")
    # %%%%%%%%%%%%%%%%
    # Create ids of studies
    ns <- vapply(seq_len(n_studies), function(study_index) {
      paste0(sample(c(letters, LETTERS, 0:9), 7L, replace = TRUE), collapse = "")
    }, character(1))
    assertthat::assert_that(length(unique(ns)) == n_studies)
    
    out <- tidyr::crossing(
      SS = paste0(ns, " ", 1),
      SSBS = 1:n_sites) |> 
      dplyr::mutate(SSBS = paste0(SS, " ", SSBS),
                    # TODO: Add blocking scheme
                    SSB = SS,
                    Source_ID = substr(SS, start = 0, stop = 7),
                    Predominant_land_use = sample(lu_legend, size = n_sites * n_studies, replace = TRUE),
                    Use_intensity = sample(ui_legend, size = n_sites * n_studies, replace = TRUE, prob = c(1,.5,.2))
                    )
    
    # Other study settings
    out$Sampling_effort <- stats::rexp(nrow(out),rate = 2)
    # Randomly sample some studies (SS) and keep effort constant (1)
    const_effort_SS <- sample(unique(out$SS), round(length(unique(out$SS))/3*2))
    out$Sampling_effort[out$SS %in% const_effort_SS] <- 1
    
    out$Diversity_metric_type <- "Abundance"
    out$Diversity_metric <- "Abundance"
    out$Diversity_metric_is_suitable_for_Chao <- TRUE
    out$Diversity_metric_is_effort_sensitive <- TRUE
    out$Higher_taxon <- "Dragon"
    
    out$Sample_start_earliest <- as.POSIXct("2024-01-01")
    out$Sample_end_latest <- as.POSIXct("2024-02-01")
    
    # Lat/Long sampling
    # get a coarse-ish world map...
    rlang::check_installed(c("sf", "rnaturalearth", "rnaturalearthdata"),
                 reason = "to generate dummy species data.")
    map <- rnaturalearth::ne_countries(scale = 110, returnclass = 'sf')
    map <- map |>
      dplyr::filter(continent != "Antarctica")
    # Issues with Russia and Sudan
    map <- sf::st_make_valid(map)
    map <- map[sf::st_is_valid(map),] # Rm Sudan
    
    ss_centres <- sf::st_sample(map, n_studies)
    ss_buffers <- sf::st_buffer(ss_centres, dist = 50000) # 50km radius - could adapt
    
    ssbs_samp <- vector("list", n_studies)
    
    for (i in 1:length(ss_buffers)){
      tmp <- NULL
      tmp_out <- NULL
      tmp <- sf::st_intersection(map, ss_buffers[i])
      tmp_out <- sf::st_coordinates(sf::st_sample(tmp, n_sites))
      ssbs_samp[[i]] <- as.data.frame(tmp_out)
    }
    
    out <- dplyr::bind_cols(out, 
                            ssbs_samp |> 
                              dplyr::bind_rows() |>
                              dplyr::rename(Latitude = Y,
                                            Longitude = X))
    # %%%%%%%%%%
    # Intersect each data point with the map (and get country, +)
    # continent, region_un, formal_en, admin

    map_samp_idx <- ssbs_samp |> 
      dplyr::bind_rows() |>
      sf::st_as_sf(coords = c(1:2), crs = sf::st_crs(map)) |>
      sf::st_intersects(map, sparse = F) |>
      apply(1, which)
      
    # Extract some geographic info, note that region_un and subregion are not 
    # quIte aligned with UN_region/UN_subregion
    extra_dat <- as.data.frame(map)[map_samp_idx,
                                    c(#"continent", 
                                      "region_un", 
                                      "subregion",  
                                      "formal_en"
                                      # "region_wb",
                                      # "admin", 
                                      # "sovereignt"
                                      ) ] |>
      dplyr::rename(UN_region = region_un,
                    UN_subregion = subregion,
                    Country = formal_en)

    out <- dplyr::bind_cols(out,
                            extra_dat)
    # %%%%%%%%%%
    
    # spp abund sampling...
    spp_ls <- vector("list", n_studies)
    cntr <- 1
    for (ss in unique(out$SS)){
      # For each SS, sample N spp - 
      ssbs <- unique(out$SSBS[out$SS == ss])
      N_spp <- stats::runif(n = 1, min = 5, max = 100) |> round()
      spp_id <- replicate(N_spp,
              paste0(sample(c(letters, LETTERS, 0:9), 7L, replace = TRUE),
                 collapse = ""))
      
      # Sample mean abundance (log) - runif()
      mean_n <- stats::runif(n = 1, min = 5, max = 100) |> round()
      abund_vec <- stats::rexp(N_spp * n_sites, rate = 1/mean_n) |> round()
      
      tmp_spp <-  tidyr::crossing(SSBS = ssbs, 
                                  Best_guess_binomial = spp_id) |>
        dplyr::mutate(Measurement = abund_vec,
                      Taxon_name_entered = Best_guess_binomial)
      
      spp_ls[[cntr]] <- tmp_spp
      cntr <- cntr + 1
    }
    
    spp_df <- spp_ls |> dplyr::bind_rows()
    
    out <- spp_df |>
      dplyr::left_join(out,
                       by = "SSBS")
    
  }
  return(out)
}
