#' Calculate measures of dissimilarity (wrapper) between sites
#' 
#' @description 
#' This function calculates dissimilarity measures between sites within a study,
#' for use as covariates in models of comopositional dis/similarity
#' @details 
#' Expects a [`HibDatabase`] object.
#' @param obj A [`HibDatabase`] object.
#' @param metrics A [`character`] vector of 'distance' metrics
#' @param baseline_lui A [`character`] or NULL, indicating the baseline LUI against
#' which compositional dis/similarity is calculated per study - should match the 
#' value used when calculating compositional dis/similarity (hib_calc_div_between)
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_site_contrasts(obj, 
#'                               c("lui_contrast", "geo_distance", "env_distance"),
#'                               T)
#' } 
#' @export
hib_calc_site_contrasts <- function(obj,
                                    metrics,
                                    baseline_lui = "Primary vegetation_Minimal use"){
  
  # TODO:
  # convert primin_baseline to baseline_lui, and lui_colname
  # lui_colname may be unnecessary, - lui_legend assumes LUI...
  
  # Assertion checks
  assertthat::assert_that(
    inherits(obj, "HibDatabase"),
    is.character(metrics),
    # is.logical(primin_baseline) 
    is.character(baseline_lui) || is.null(baseline_lui)
  )
  
  # check metrics are valid options
  metrics <- match.arg(metrics,
                       c("lui_contrast", 
                         "geogr_dist_hav", "geogr_dist_geod",
                         "geogr_midpoint"),
                       # , "env_distance"
                       several.ok = T)
  
  message(crayon::green("Collating data for site contrasts...\n"))
  
  # get obj$sites columns of interest - SSBS, SS, LUI, lat/long etc
  cols_of_interest <- c("SSBS", "SSB", "SS")
  if ("lui_contrast" %in% metrics){
    cols_of_interest <- c(cols_of_interest, "LUI")
  }
  if (any(c("geogr_dist_hav", "geogr_dist_geod", "geogr_midpoint") %in% metrics)){
    cols_of_interest <- c(cols_of_interest, "Latitude", "Longitude")
  }
  
  
  # assert that cols of interest are in obj$sites
  # subset of c("SSBS", "SSB", "SS", "LUI", "Latitude", "Longitude")
  assertthat::assert_that(
    all(cols_of_interest %in% colnames(obj$sites))
  )
  site_df <- obj$sites[, cols_of_interest]
  
  # %%%
  # Get site contrasts

  # assert that specific cols are in obj$species 
  # c(Predominant_land_use, Diversity_metric_type, SS, Sampling_effort, Taxon_name_entered)
  
  assertthat::assert_that(
    all(c("Predominant_land_use", "LUI", "Diversity_metric_type", "SS", 
          "SSBS", "Sampling_effort", "Taxon_name_entered") %in% colnames(obj$species))
  )
  # NOTE: they should be anyway based on hib database self checks...

  spp_tmp <- obj$species
  
  spp_tmp <- spp_tmp |>
    dplyr::filter(!is.na(Predominant_land_use) &
                    Diversity_metric_type == "Abundance") |>
    dplyr::group_by(SS) |>
    dplyr::mutate(n_sample_effort = dplyr::n_distinct(Sampling_effort),
           # Taxon_name_entered v best_guess_binomial ??
           n_species = dplyr::n_distinct(Taxon_name_entered)) |>
    dplyr::ungroup() |>
    dplyr::filter(n_sample_effort == 1 &
             n_species > 1)
  
  # if (primin_baseline == T){
  #   # TODO: should this be ammended to primary forest...
  #   spp_tmp$prim_min <- as.numeric((grepl("primary", spp_tmp$Predominant_land_use, ignore.case = T) &
  #                                     grepl("minimal", spp_tmp$Use_intensity, ignore.case = T)))
  #   spp_tmp <- spp_tmp |>
  #     dplyr::group_by(SS) |>
  #     dplyr::mutate(n_primin_records = sum(prim_min)) |>
  #     dplyr::filter(n_primin_records > 0) |>
  #     dplyr::ungroup()
  #   
  #   site_tmp <- obj$sites |>
  #     dplyr::filter(SSBS %in% unique(spp_tmp$SSBS)) |>
  #     left_join(unique(spp_tmp[, c("SSBS", "prim_min")])) |>
  #     droplevels()
  #   
  # } else if (primin_baseline == F){
  #   site_tmp <- obj$sites |>
  #     dplyr::filter(SSBS %in% unique(spp_tmp$SSBS)) |>
  #     droplevels()
  # }
  if (!is.null(baseline_lui)){
    assertthat::assert_that(baseline_lui %in% unique(spp_tmp$LUI))

    spp_tmp$base_lui <- spp_tmp$LUI == baseline_lui
    
    spp_tmp <- spp_tmp |>
      dplyr::group_by(SS) |>
      dplyr::mutate(n_base_lui = sum(base_lui)) |>
      dplyr::filter(n_base_lui > 0) |>
      dplyr::ungroup()
    
    site_tmp <- obj$sites |>
      dplyr::filter(SSBS %in% unique(spp_tmp$SSBS)) |>
      dplyr::left_join(unique(spp_tmp[, c("SSBS", "base_lui")])) |>
      droplevels()
    
  } else if (is.null(baseline_lui)){
    site_tmp <- obj$sites |>
      dplyr::filter(SSBS %in% unique(spp_tmp$SSBS)) |>
      droplevels()
  }
  
  # Need to get combinations of data... - s1, s2, ...
  # %%%
  # From https://adrianadepalma.github.io/BII_tutorial/bii_example.html
  
  # get a vector of each study to loop over
  studies <- unique(as.character(site_tmp$SS))
  
  site_comp_df <- purrr::map_dfr(.x = studies, .f = function(x){
    
    # pull out the sites that are Primary minimal (we only want to use comparisons with the baseline)
    # if (primin_baseline == T){
    #   site_data <- dplyr::filter(site_tmp, SS == x) |>
    #     # pull out the SSBS and LandUse information
    #     dplyr::select(SSBS, prim_min)
    #   
    #   baseline_sites <- site_data |>
    #     dplyr::filter(prim_min == 1) |>
    #     dplyr::pull(SSBS)
    # } else if (primin_baseline == F){
    #   site_data <- dplyr::filter(site_tmp, SS == x)
    #   
    #   baseline_sites <- site_data |>
    #     dplyr::pull(SSBS)
    # }
    if (!is.null(baseline_lui)){
      site_data <- dplyr::filter(site_tmp, SS == x) |>
        # pull out the SSBS and LandUse information
        dplyr::select(SSBS, base_lui)
      
      baseline_sites <- site_data |>
        dplyr::filter(base_lui == T) |>
        dplyr::pull(SSBS)
    } else if (is.null(baseline_lui)){
      site_data <- dplyr::filter(site_tmp, SS == x)
      
      baseline_sites <- site_data |>
        dplyr::pull(SSBS)
    }
    
    # pull out all the sites
    site_list <- site_data |>
      dplyr::pull(SSBS)
    
    # get all site x site comparisons for this study
    site_comparisons <- expand.grid(baseline_sites, site_list) |>
      # rename the columns so they will be what the compositional similarity function expects for ease
      dplyr::rename(SSBS1 = Var1, SSBS2 = Var2) |>
      # make the values characters rather than factors
      dplyr::mutate(SSBS1 = as.character(SSBS1),
                    SSBS2 = as.character(SSBS2)) |>
      # remove the comparisons where the same site is being compared to itself
      dplyr::filter(SSBS1 != SSBS2) |>
      # add the study id
      dplyr::mutate(SS = as.character(x)) 
    
    return(site_comparisons)
  })
  
  message(crayon::green(paste0(nrow(site_comp_df), 
                               " site comparisons to make...\n"))
  )
  
  message(crayon::green("Calculating site contrasts...\n"))
  
  # %%%
  # calculate metrics of interest
  # - LUI - left join SSBS1, LUI1, SSBS2, LUI2, paste together
  if ("lui_contrast" %in% metrics){
    # Edited from https://adrianadepalma.github.io/BII_tutorial/bii_example.html
    site_comp_df <- site_comp_df |>
      dplyr::left_join(site_df[, c("SSBS", "LUI")], 
                       by = c("SSBS1" = "SSBS")) |>
      dplyr::rename(SSBS1_LUI = LUI) |>
      # get the land use for s2
      dplyr::left_join(site_df[, c("SSBS", "LUI")], 
                       by = c("SSBS2" = "SSBS")) |>
      dplyr::rename(SSBS2_LUI = LUI) |>
      # create an lui_contrast column (what we'll use for modelling)
      dplyr::mutate(LUI_contrast = paste(SSBS1_LUI, SSBS2_LUI, sep = "_vs_"),
                    Baseline_LUI = baseline_lui)
    if (is.null(baseline_lui)){
      site_comp_df <- site_comp_df |>
        dplyr::mutate(Baseline_LUI = NA)
    }
  }
  # TODO: lu contrast (Predominant_land_use), i contrast (Use_intensity)
  # if ("lu_contrast" %in% metrics){
  #   # Edited from https://adrianadepalma.github.io/BII_tutorial/bii_example.html
  #   site_comp_df <- site_comp_df |>
  #     dplyr::left_join(site_df[, c("SSBS", "Predominant_land_use")], 
  #                      by = c("SSBS1" = "SSBS")) |>
  #     dplyr::rename(SSBS1_LU = Predominant_land_use) |>
  #     # get the land use for s2
  #     dplyr::left_join(site_df[, c("SSBS", "Predominant_land_use")], 
  #                      by = c("SSBS2" = "SSBS")) |>
  #     dplyr::rename(SSBS2_LU = Predominant_land_use) |>
  #     # create an lu_contrast column (what we'll use for modelling)
  #     dplyr::mutate(LU_contrast = paste(SSBS1_LU, SSBS2_LU, sep = "_vs_"),
  #                   Baseline_LUI = baseline_lui) # ??
  #   if (is.null(baseline_lui)){  # ?? Baseline_LUI isn't really needed for LU??
  #     site_comp_df <- site_comp_df |>
  #       dplyr::mutate(Baseline_LUI = NA)
  #   }
  # }
  # if ("i_contrast" %in% metrics){
  #   # Edited from https://adrianadepalma.github.io/BII_tutorial/bii_example.html
  #   site_comp_df <- site_comp_df |>
  #     dplyr::left_join(site_df[, c("SSBS", "Use_intensity")], 
  #                      by = c("SSBS1" = "SSBS")) |>
  #     dplyr::rename(SSBS1_I = Use_intensity) |>
  #     # get the land use for s2
  #     dplyr::left_join(site_df[, c("SSBS", "Use_intensity")], 
  #                      by = c("SSBS2" = "SSBS")) |>
  #     dplyr::rename(SSBS2_I = Use_intensity) |>
  #     # create an lu_contrast column (what we'll use for modelling)
  #     dplyr::mutate(I_contrast = paste(SSBS1_I, SSBS2_I, sep = "_vs_"),
  #                   Baseline_LUI = baseline_lui)
  #   if (is.null(baseline_lui)){
  #     site_comp_df <- site_comp_df |>
  #       dplyr::mutate(Baseline_LUI = NA)
  #   }
  # }

  # - geo_distance - SSBS1, long1, lat1, SSBS2, long2, lat2
  if (any(c("geogr_dist_hav", "geogr_dist_geod", "geogr_midpoint") %in% metrics)){
    # Edited from https://adrianadepalma.github.io/BII_tutorial/bii_example.html
    site_comp_df <- site_comp_df |>
      # get the lat and long for s1
      dplyr::left_join(site_df[, c("SSBS", "Longitude" ,"Latitude")], 
                       by = c("SSBS1" = "SSBS")) |>
      dplyr::rename(SSBS1_latitude = Latitude,
                    SSBS1_longitude = Longitude) |>
      # get the lat and long for s2
      dplyr::left_join(site_df[, c("SSBS", "Longitude" ,"Latitude")], 
                       by = c("SSBS2" = "SSBS")) |>
      dplyr::rename(SSBS2_latitude = Latitude,
                    SSBS2_longitude = Longitude)
    
    if ("geogr_dist_hav" %in% metrics){
      site_comp_df <- site_comp_df |>
        # calculate the geographic distances between s1 and s2 
        dplyr::mutate(Geogr_dist_hav = geosphere::distHaversine(cbind(SSBS1_longitude, SSBS1_latitude), 
                                                                cbind(SSBS2_longitude, SSBS2_latitude)))
    }
    if ("geogr_dist_geod" %in% metrics){
      site_comp_df <- site_comp_df |>
        # calculate the geographic distances between s1 and s2 
        dplyr::mutate(Geogr_dist_geod = geosphere::distGeo(cbind(SSBS1_longitude, SSBS1_latitude), 
                                                           cbind(SSBS2_longitude, SSBS2_latitude)))
    }
    if ("geogr_midpoint" %in% metrics){
      site_comp_df <- site_comp_df |>
        # calculate the geographic midpoint between s1 and s2 
        dplyr::mutate(Latitude_midpoint = (SSBS1_latitude + SSBS2_latitude) / 2,
                      Longitude_midpoint = (SSBS1_longitude + SSBS2_longitude) / 2)
    }
      
  }
  # Env dist? - gower dist [altitude, bioclim: maximum and minimum temperature and precipitation of wettest and driest month]
  
  # Also, by default incl SSB - either SSB or NA if sites are from different blocks
  site_comp_df <- site_comp_df |>
    # get the lat and long for s1
    dplyr::left_join(site_df[, c("SSBS", "SSB")], 
                     by = c("SSBS1" = "SSBS")) |>
    dplyr::rename(SSB1 = SSB) |>
    dplyr::left_join(site_df[, c("SSBS", "SSB")], 
                     by = c("SSBS2" = "SSBS")) |>
    dplyr::rename(SSB2 = SSB) |>
    dplyr::mutate(SSB = dplyr::case_when(SSB1 == SSB2 ~ SSB1,
                                         T ~ NA))
  
  # TODO:
  # add in version info... ?
  
  # %%%
  # add site_contrasts to obj
  site_contr_empty <- is.null(obj$site_contrasts)
  
  if (!site_contr_empty){
    obj$site_contrasts <- dplyr::left_join(obj$site_contrasts, 
                                           site_comp_df)
  }  
  if (site_contr_empty){
    obj$site_contrasts <- site_comp_df
  }
  
  obj$last_update <- lubridate::now()
  
  return(obj)

}
