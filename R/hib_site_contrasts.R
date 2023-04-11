#' Calculate measures of dissimilarity (wrapper) between sites
#' 
#' @description 
#' This function calculates dissimilarity measures between sites within a study,
#' for use as covariates in models of comopositional dis/similarity
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @param metrics A [`character`] vector of dissimilarity metrics...
#' @param primin_baseline A [`logical`] value, determines whether 
#' dissimilarity is calculated for all site combinations per study, or
#' only relative to a primary minimal baseline - should match the value used when
#' caclulating compositional dis/similarity (hib_calc_div_between)
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
                                    primin_baseline = T){
  
  # TODO:
  # convert primin_baseline to baseline_lui
  # Assertion checks
  assertthat::assert_that(
    is.character(metrics),
    is.logical(primin_baseline) 
  )
  
  # check metrics are valid options
  metrics <- match.arg(metrics,
                       c("lui_contrast", 
                         "geogr_dist_hav", "geogr_dist_geod"),
                       # , "env_distance"
                       several.ok = T)
  
  message(crayon::green("Collating data for site contrasts...\n"))
  
  # get obj$sites columns of interest - SSBS, SS, LUI, lat/long etc
  cols_of_interest <- c("SSBS", "SSB", "SS")
  if ("lui_contrast" %in% metrics){
    cols_of_interest <- c(cols_of_interest, "LUI")
  }
  if (any(c("geogr_dist_hav", "geogr_dist_geod") %in% metrics)){
    cols_of_interest <- c(cols_of_interest, "Latitude", "Longitude")
  }
  
  site_df <- obj$sites[, cols_of_interest]
  
  # %%%
  # Get site contrasts
  spp_tmp <- obj$species
  
  spp_tmp <- spp_tmp %>%
    dplyr::filter(!is.na(Predominant_land_use) &
                    Diversity_metric_type == "Abundance") %>%
    group_by(SS) %>%
    mutate(n_sample_effort = n_distinct(Sampling_effort),
           # Taxon_name_entered v best_guess_binomial ??
           n_species = n_distinct(Taxon_name_entered)) %>%
    ungroup() %>%
    filter(n_sample_effort == 1 &
             n_species > 1)
  
  if (primin_baseline == T){
    # TODO: should this be ammended to primary forest...
    spp_tmp$prim_min <- as.numeric((grepl("primary", spp_tmp$Predominant_land_use, ignore.case = T) &
                                      grepl("minimal", spp_tmp$Use_intensity, ignore.case = T)))
    spp_tmp <- spp_tmp %>%
      dplyr::group_by(SS) %>%
      dplyr::mutate(n_primin_records = sum(prim_min)) %>%
      dplyr::filter(n_primin_records > 0) %>%
      dplyr::ungroup()
    
    site_tmp <- obj$sites %>%
      dplyr::filter(SSBS %in% unique(spp_tmp$SSBS)) %>%
      left_join(unique(spp_tmp[, c("SSBS", "prim_min")])) %>%
      droplevels()
    
  } else if (primin_baseline == F){
    site_tmp <- obj$sites %>%
      dplyr::filter(SSBS %in% unique(spp_tmp$SSBS)) %>%
      droplevels()
  }
  
  
  # n_spp <- spp_tmp %>%
  #   dplyr::select(Best_guess_binomial) %>%
  #   dplyr::n_distinct()
  
  # Need to get combinations of data... - s1, s2, ...
  # %%%
  # From https://adrianadepalma.github.io/BII_tutorial/bii_example.html
  
  # get a vector of each study to loop over
  studies <- unique(as.character(site_tmp$SS))
  
  site_comp_df <- purrr::map_dfr(.x = studies, .f = function(x){
    
    # pull out the sites that are Primary minimal (we only want to use comparisons with the baseline)
    if (primin_baseline == T){
      site_data <- dplyr::filter(site_tmp, SS == x) %>%
        # pull out the SSBS and LandUse information
        dplyr::select(SSBS, prim_min)
      
      baseline_sites <- site_data %>%
        dplyr::filter(prim_min == 1) %>%
        dplyr::pull(SSBS)
    } else if (primin_baseline == F){
      site_data <- dplyr::filter(site_tmp, SS == x)
      
      baseline_sites <- site_data %>%
        dplyr::pull(SSBS)
    }
    
    # pull out all the sites
    site_list <- site_data %>%
      dplyr::pull(SSBS)
    
    # get all site x site comparisons for this study
    site_comparisons <- expand.grid(baseline_sites, site_list) %>%
      # rename the columns so they will be what the compositional similarity function expects for ease
      dplyr::rename(SSBS1 = Var1, SSBS2 = Var2) %>%
      # remove the comparisons where the same site is being compared to itself
      dplyr::filter(SSBS1 != SSBS2) %>%
      # make the values characters rather than factors
      dplyr::mutate(SSBS1 = as.character(SSBS1),
                    SSBS2 = as.character(SSBS2),
                    # add the full name
                    # contrast = paste(SSBS1, SSBS2, sep = "_"),
                    # add the study id
                    SS = as.character(x)) 
    
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
    site_comp_df <- site_comp_df %>%
      dplyr::left_join(site_df[, c("SSBS", "LUI")], 
                       by = c("SSBS1" = "SSBS")) %>%
      dplyr::rename(SSBS1_LUI = LUI) %>%
      # get the land use for s2
      dplyr::left_join(site_df[, c("SSBS", "LUI")], 
                       by = c("SSBS2" = "SSBS")) %>%
      dplyr::rename(SSBS2_LUI = LUI) %>%
      # create an lu_contrast column (what we'll use for modelling)
      dplyr::mutate(LUI_contrast = paste(SSBS1_LUI, SSBS2_LUI, sep = "_vs_"))
  }
  # - geo_distance - SSBS1, long1, lat1, SSBS2, long2, lat2
  if (any(c("geogr_dist_hav", "geogr_dist_geod") %in% metrics)){
    # Edited from https://adrianadepalma.github.io/BII_tutorial/bii_example.html
    site_comp_df <- site_comp_df %>%
      # get the lat and long for s1
      dplyr::left_join(site_df[, c("SSBS", "Longitude" ,"Latitude")], 
                       by = c("SSBS1" = "SSBS")) %>%
      dplyr::rename(SSBS1_latitude = Latitude,
                    SSBS1_longitude = Longitude) %>%
      # get the lat and long for s2
      dplyr::left_join(site_df[, c("SSBS", "Longitude" ,"Latitude")], 
                       by = c("SSBS2" = "SSBS")) %>%
      dplyr::rename(SSBS2_latitude = Latitude,
                    SSBS2_longitude = Longitude)
    
    if ("geogr_dist_hav" %in% metrics){
      site_comp_df <- site_comp_df %>%
        # calculate the geographic distances between s1 and s2 sites
        dplyr::mutate(geogr_dist_hav = geosphere::distHaversine(cbind(SSBS1_longitude, SSBS1_latitude), 
                                                                cbind(SSBS2_longitude, SSBS2_latitude)))
    }
    if ("geogr_dist_geod" %in% metrics){
      site_comp_df <- site_comp_df %>%
        # calculate the geographic distances between s1 and s2 sites
        dplyr::mutate(geogr_dist_geod = geosphere::distGeo(cbind(SSBS1_longitude, SSBS1_latitude), 
                                                           cbind(SSBS2_longitude, SSBS2_latitude)))
    }
      
  }
  
  # Also, by default incl SSB - either SSB or NA if sites are from different blocks
  site_comp_df <- site_comp_df %>%
    # get the lat and long for s1
    dplyr::left_join(site_df[, c("SSBS", "SSB")], 
                     by = c("SSBS1" = "SSBS")) %>%
    dplyr::rename(SSB1 = SSB) %>%
    dplyr::left_join(site_df[, c("SSBS", "SSB")], 
                     by = c("SSBS2" = "SSBS")) %>%
    dplyr::rename(SSB2 = SSB) %>%
    dplyr::mutate(SSB = dplyr::case_when(SSB1 == SSB2 ~ SSB1,
                                         T ~ NA))
  
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
  
  return(obj)

}
