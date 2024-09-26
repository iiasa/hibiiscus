#' Calculate measures of dissimilarity (wrapper) between sites
#' 
#' @description 
#' This function calculates dissimilarity measures between sites within a study,
#' for use as covariates in models of comopositional dis/similarity
#' @details 
#' Expects a [`HibDatabase`] object.
#' @param obj A [`HibDatabase`] object.
#' @param metrics A [`character`] vector of 'distance' metrics. Current options are;
#' "lui_contrast", "geogr_dist_hav", "geogr_dist_geod" and "geogr_midpoint", "env_dist".
#' @param custom_metric_specs A named [`list`] with each named item containing the
#' specifications for the named custom metric. Custom metrics must be based on 
#' vegdist distance metrics. 
#' Essential named values within the specifications of a custom metric are: 
#' \code{dist_metr} - a character indicating the metric of choice;
#' \code{vars} - a character vector indicating the site variables to use when calculating
#' the distance metric. 
#' For example: \code{list(Gower_dist_a = list(dist_metr = "gower", vars = c("env1", "env2", "env3")),
#'                         Gower_dist_b = list(dist_metr = "gower", vars = c("env4", "env5", "env6")))}
#' Metric names (as specified in the name of the sub-lists) will be converted to 
#' sentence-case in the output data.
#' @param baseline_lui A [`character`] or NULL, indicating the baseline LUI against
#' which dis/similarity is calculated per study - should match the value used 
#' when calculating compositional dis/similarity (hib_calc_div_between)
#' @param overwrite A [`logical`] indicating whether or not any current site contrast
#' data should be over-written.
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_site_contrasts(obj, 
#'                               c("lui_contrast", "geogr_dist_hav"),
#'                               "Primary vegetation_Minimal use")
#' } 
#' @export
hib_calc_site_contrasts <- function(obj,
                                    metrics,
                                    custom_metric_specs = NULL, #list(list(custom_name = "...", dist_metr = "gower", vars = c("...")), list())
                                    baseline_lui = "Primary vegetation_Minimal use",
                                    overwrite = FALSE){
  
  # TODO:
  # lui_colname may be unnecessary, - lui_legend assumes LUI...
  
  # TODO:
  # Custom metric - 
  # use "custom_metric" in metrics vec
  # add a custom_metric_spec arg - list/named vector...
  # have custom_metric_func, custom_metric_max, custom_metric_name 
  # How can this work...
  
  # Assertion checks
  assertthat::assert_that(
    inherits(obj, "HibDatabase"),
    is.character(metrics),
    is.character(baseline_lui) || is.null(baseline_lui),
    is.list(custom_metric_specs) || is.null(custom_metric_specs)
  )
  
  # check metrics are valid options
  metrics <- match.arg(metrics,
                       c("lui_contrast", 
                         "geogr_dist_hav", "geogr_dist_geod",
                         "geogr_midpoint", "env_dist",
                         "custom_dist"),
                       # 
                       # "custom_dist"
                       several.ok = T)
  # TODO:
  # Custom dist specs... 
  # make a named list - i.e., list(Custom_dist1 = list(...), Custom_dist_2 = list(...)))
  # Much neater!!
  # Name, dist_type (vegan), covar cols

  if (overwrite){
    obj$site_contrasts <- NULL
  }

  
  # get obj$sites columns of interest - SSBS, SS, LUI, lat/long etc
  cols_of_interest <- c("SSBS", "SSB", "SS")
  if ("lui_contrast" %in% metrics){
    cols_of_interest <- c(cols_of_interest, "LUI")
  }
  if (any(c("geogr_dist_hav", "geogr_dist_geod", "geogr_midpoint") %in% metrics)){
    cols_of_interest <- c(cols_of_interest, "Latitude", "Longitude")
  }
  if ("env_dist" %in% metrics){
    cols_of_interest <- c(cols_of_interest, "wc2.1_30s_bio_5", "wc2.1_30s_bio_6", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", "wc2.1_30s_elev")
  }

  # %%%%%%
  # TODO:
  # Custom distance calcs...
  if ("custom_dist" %in% metrics){  
    assertthat::assert_that(
      is.list(custom_metric_specs) 
      )

    # TODO:
    # assert that all items of list have correct names /formats...
    assertthat::assert_that(
      lapply(custom_metric_specs,
             function(x){ all( c("dist_metr", "vars") %in% names(x) )}) |>
      unlist() |>
      all()
      )
    
    # custom_name , dist_metr, vars

    # all present, all chr

    # Get custom names
    custom_metric_names <- names(custom_metric_specs)
    # lapply(custom_metric_specs, function(x){x$custom_name}) |> unlist()

    default_metric_names <- c("LUI_contrast", 
                              "Geogr_dist_hav", "Geogr_dist_geod",
                              "Geogr_midpoint", "Env_dist")
    # Check names of metrics (no-overlap with defaults.. - if so, add "_custom", and vars of interest 
    if (any(custom_metric_names %in% default_metric_names)){
      # Alternatively, just stop?
      # TODO: Warning message
      edit_idx <- which(custom_metric_names %in% default_metric_names)

      old_nms <- names(custom_metric_specs)[edit_idx]

      names(custom_metric_specs)[edit_idx] <- paste0(names(custom_metric_specs)[edit_idx], "_custom")
      # for (idx in edit_idx){
      #   custom_metric_specs[[idx]]$custom_name <- paste0(custom_metric_specs[[idx]]$custom_name, "_custom")
      # }
      new_nms <- names(custom_metric_specs)[edit_idx]

      prnt_df <- data.frame(old_nms = old_nms,
                            new_nms = new_nms) |>
        dplyr::mutate(prnt_out = paste(old_nms, new_nms, sep = " -->> "))

      message(crayon::yellow(paste0("\nWarning: custom metric names match default ",
                                  "contrasts! The following names have been ",
                                  "adjusted;\n\t",
                                  paste(prnt_df$prnt_out, collapse = "\n\t")
                                  # "\n"
                                  )))
    }

    custom_metric_names <- names(custom_metric_specs)
    # lapply(custom_metric_specs, function(x){x$custom_name}) |> unlist()
    
    # print(custom_metric_names)

    if (length(unique(custom_metric_names)) < length(custom_metric_names)){
      # print("duplicates")
      # Alternatively, just stop?
      old_nms <- custom_metric_names
      dupl_nms <- custom_metric_names[duplicated(custom_metric_names)]

      # TODO: Warning message
      # For unique names, id N, add _1:n to end of name...
      un_dupl_nms <- unique(dupl_nms)
      for (nm in un_dupl_nms){
        reps <- which(custom_metric_names == nm)
        n_reps <- length(reps)

        for (i in 1:n_reps){

          names(custom_metric_specs)[reps[i]] <- paste0(names(custom_metric_specs)[reps[i]], "_", as.character(i))
          # custom_metric_specs[[reps[i]]]$custom_name <- paste0(custom_metric_specs[[reps[i]]]$custom_name, "_", as.character(i))
        }
      }
      new_nms <- names(custom_metric_specs)

      # %%%%%
      prnt_df <- data.frame(old_nms = old_nms,
                            new_nms = new_nms) |>
        dplyr::filter(old_nms != new_nms) |>
        dplyr::mutate(prnt_out = paste(old_nms, new_nms, sep = " -->> "))

      message(crayon::yellow(paste0("\nWarning: multiple custom metric names are the same! ",
                                  "The following names have been adjusted;\n\t",
                                  paste(prnt_df$prnt_out, collapse = "\n\t")
                                  # "\n"
                                  )))
      # %%%%%
      # print(custom_metric_specs)

      # Finally, what about if any custom_metric_names are already in the site_contrasts df - ignore?
      # Stop?

    }


    # custom_metric_dist <- lapply(custom_metric_specs, function(x){x$dist_metr}) |> unlist() |> unique()
    
    # Get custom vars
    custom_metric_vars <- lapply(custom_metric_specs, function(x){x$vars}) |> unlist() |> unique()

    cols_of_interest <- c(cols_of_interest, custom_metric_vars)
  }

  
  # assert that cols of interest are in obj$sites
  # subset of c("SSBS", "SSB", "SS", "LUI", "Latitude", "Longitude")
  assertthat::assert_that(
    all(cols_of_interest %in% colnames(obj$sites))
  )
  site_df <- obj$sites[, cols_of_interest]
  
  # %%%
  # Get site contrasts
  message(crayon::green("\nCollating data for site contrasts...\n"))
  
  # TODO: if site_contrasts already exists in hib object, use this [SS, SSBS...]
  
  # option to incl over-write...
  # site_contr_empty <- is.null(obj$site_contrasts)
  if ( (!is.null(obj$site_contrasts)) & (overwrite == F) ){
    site_comp_df <- obj$site_contrasts[, c("SSBS1", "SSBS2", "SS")]
    message(crayon::green(paste0("\t", nrow(site_comp_df), 
                                 " site comparisons to make...\n"))
    )
    # Check if metrics are already present! - if they are, return a warning/error,
    metric_out_names <- tools::toTitleCase(metrics)
    metric_out_names[metric_out_names == "Lui_contrast"] <- "LUI_contrast"
    metric_out_names <- metric_out_names[metric_out_names != "Custom_dist"]
    
    new_metrics_idx <- which(!(metric_out_names %in% colnames(obj$site_contrasts)))
    # drop from metrics and print the dropped metrics...
    metrics <- metrics[new_metrics_idx]
    dropped_metrics <- metric_out_names[metric_out_names %in% colnames(obj$site_contrasts)]
    

    # Also check custom metric names (names(custom_metric_specs))
    if (!is.null(custom_metric_specs)){
      if (any(names(custom_metric_specs) %in% colnames(obj$site_contrasts))){ 
        # drop them!, add to dropped metrics
        keep_idx <- which(!names(custom_metric_specs) %in% colnames(obj$site_contrasts))
        drop_idx <- which(names(custom_metric_specs) %in% colnames(obj$site_contrasts))
        dropped_custom_metrics <- names(custom_metric_specs)[drop_idx]
        dropped_metrics <- c(dropped_metrics, dropped_custom_metrics)

        custom_metric_specs <- custom_metric_specs[[keep_idx]]
        if (length(custom_metric_specs) == 0){
          custom_metric_specs <- NULL
        }
      }
    }

    if (length(dropped_metrics) > 0){
      message(crayon::yellow(paste0("\nWarning: the following contrasts have already been ",
                                  "added to the HibDatabase and so will not be ",
                                  "recalculated;\n\t",
                                  paste(dropped_metrics, collapse = ", ")
                                  # "\n"
                                  )))
    }
    # if nrow(metrics_spec_df) == 0
    # print(nrow(metrics_spec_df))
    # TODO: what about custom metrics...!!!
    metrics_prnt <- tools::toTitleCase(metrics)
    metrics_prnt[metrics_prnt == "Lui_contrast"] <- "LUI_contrast"
    metrics_prnt <- metrics_prnt[metrics_prnt != "Custom_dist"]

    if (!is.null(custom_metric_specs)){
      metrics_prnt <- c(metrics_prnt, names(custom_metric_specs))
    }
    # use metrics_prnt

    if (length(metrics_prnt) == 0){
        # message(crayon::red(paste0("\nError: no uncalculated biodiversity metrics.\n\n",
        #                        "Check that the specified biodiversity responses ('", baseline_lui,
        #                        "') is present ...\n"))
        #   )
        stop("No new contrast metrics to calculate.", call. = F)
    }

    # TODO Check for same baseline LUI??

  } 

  if ( (is.null(obj$site_contrasts)) | (overwrite == T) ){
    # If site contrasts doesn't exist/over-write is True - generate site_comp_df from scratch
    # Set contrasts to NULL?? - should over-write over-write the whole site_contrasts df, or just the metrics of interest... - prob the whole thing
    # %%%%%%%%%
    # assert that specific cols are in obj$species 
    # c(Predominant_land_use, Diversity_metric_type, SS, Sampling_effort, Taxon_name_entered)
    
    assertthat::assert_that(
      all(c("Diversity_metric_type", "SS", 
            "SSBS", "Sampling_effort", "Taxon_name_entered") %in% colnames(obj$species)),
      assertthat::has_name(obj$sites, "SSBS"),
      assertthat::has_name(obj$sites, "SS"),
      assertthat::has_name(obj$sites, "LUI")
    )
    # NOTE: they should be anyway based on hib database self checks...

    spp_tmp <- obj$species
    
    spp_tmp <- spp_tmp |>
      dplyr::filter(#!is.na(Predominant_land_use) &
                    Diversity_metric_type == "Abundance") |>
      dplyr::group_by(SS) |>
      dplyr::mutate(n_sample_effort = dplyr::n_distinct(Sampling_effort),
                    # Taxon_name_entered v best_guess_binomial ?? - what if none..
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

      # Better appr is to first get site data..., (subset to spp_tmp_SSBS) 
      # then check for baseline LUI...

      # %%%%
      site_tmp <- obj$sites |>
        dplyr::filter(SSBS %in% unique(spp_tmp$SSBS),
                    !is.na(LUI))

      # print(baseline_lui)
      # print(unique(site_tmp$LUI))
      assertthat::assert_that(baseline_lui %in% unique(site_tmp$LUI))

      site_tmp$base_lui <- site_tmp$LUI == baseline_lui
      
      site_tmp <- site_tmp |>
        dplyr::group_by(SS) |>
        dplyr::mutate(n_base_lui = sum(base_lui)) |>
        dplyr::filter(n_base_lui > 0) |>
        dplyr::ungroup() |>
        droplevels()
      # %%%%
      
      # This is an issue - what if LUI updates have only been added to the site, not spp data!!

      # assertthat::assert_that(baseline_lui %in% unique(spp_tmp$LUI))

      # spp_tmp$base_lui <- spp_tmp$LUI == baseline_lui
      
      # # TODO: 

      # spp_tmp <- spp_tmp |>
      #   dplyr::group_by(SS) |>
      #   dplyr::mutate(n_base_lui = sum(base_lui)) |>
      #   dplyr::filter(n_base_lui > 0) |>
      #   dplyr::ungroup()
      
      # site_tmp <- obj$sites |>
      #   dplyr::filter(SSBS %in% unique(spp_tmp$SSBS)) |>
      #   dplyr::left_join(unique(spp_tmp[, c("SSBS", "base_lui")]),
      #                     by = "SSBS") |>
      #   droplevels()
      
    } else if (is.null(baseline_lui)){
      site_tmp <- obj$sites |>
        dplyr::filter(SSBS %in% unique(spp_tmp$SSBS),
                    !is.na(LUI)) |>
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

    if (nrow(site_comp_df) == 0){
      message(crayon::red(paste0("\nError: no site comparisons possible...\n\n",
                                 "Maybe check that the specified baseline LUI ('", baseline_lui,
                                 "') is present ...\n"))
      )
      stop("Unable to generate any site comparisons, see above for details.")
    }
    message(crayon::green(paste0("\t", nrow(site_comp_df), 
                                 " site comparisons to make...\n"))
    )


    metrics_prnt <- tools::toTitleCase(metrics)
    metrics_prnt[metrics_prnt == "Lui_contrast"] <- "LUI_contrast"
    metrics_prnt <- metrics_prnt[metrics_prnt != "Custom_dist"]

    if (!is.null(custom_metric_specs)){
      metrics_prnt <- c(metrics_prnt, names(custom_metric_specs))
    }
  }
  

  
  # %%%%%%%%%
  
  
  
  message(crayon::green("Calculating site contrasts...\n"))
  
  # TODO: Get names in title case, with custom metric names...
  message(crayon::white(paste0("\t(",
                                 paste(metrics_prnt, 
                                       collapse = ", "),
                                 ") ...\n")))
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
      # site_comp_df <- site_comp_df |>
      #   # calculate the geographic midpoint between s1 and s2 
      #   # TODO: is therea better way to do this? - if coords are not 'true' lat/long...
      #   # then this may be a poor estimate... 
      #   # If coords are far apart, this is bad, if coords are either side of 180, this is wrong..
      #   # 
      #   # sf::st_centroid
      #   dplyr::mutate(Latitude_midpoint = (SSBS1_latitude + SSBS2_latitude) / 2,
      #                 Longitude_midpoint = (SSBS1_longitude + SSBS2_longitude) / 2)

      midpoint_df <- geosphere::midPoint(site_comp_df[, c("SSBS1_longitude", "SSBS1_latitude")], 
                                         site_comp_df[, c("SSBS2_longitude", "SSBS2_latitude")]) |>
        as.data.frame()

      colnames(midpoint_df) <- c("Longitude_midpoint", "Latitude_midpoint")
      site_comp_df <- cbind(site_comp_df, midpoint_df)
    }
      
  }



  if ("env_dist" %in% metrics){
    ss_un <- unique(site_comp_df$SS)

    site_comp_df <- site_comp_df |>
      dplyr::mutate(SSBS_contrast = paste(SSBS1, SSBS2, sep = "_"))

    # Env dist? - gower dist [altitude, bioclim: maximum and minimum temperature and precipitation of wettest and driest month]
    env_dist_df <- purrr::map_dfr(.x = ss_un,
                                 .f = function(ss = .x,
                                               env_data = site_df[,c("SS", "SSBS", "wc2.1_30s_bio_5", "wc2.1_30s_bio_6", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", "wc2.1_30s_elev")] |> 
                                                 stats::na.omit(),
                                               contrasts = site_comp_df$SSBS_contrast){
                                   
                                   tmp_env <- env_data |> 
                                     dplyr::filter(SS == ss) |> 
                                     dplyr::select(-SS) |> 
                                     tibble::column_to_rownames("SSBS")
                                   
                                   # print(tmp_env)
                                   # Need to check that this would work even if
                                   # no env data for a whole ss.. 
                                   if (nrow(tmp_env) > 1){
                                     tmp_dist <- hib_calc_veg(tmp_env,
                                                              "gower",
                                                              F,
                                                              F,
                                                              F)
                                     
                                     tmp_dist <- tmp_dist |>
                                       dplyr::rename(SSBS_contrast = contrast) |>
                                       dplyr::filter(SSBS_contrast %in% contrasts)
                                     tmp_dist$SS <- ss
                                   } else {tmp_dist <- NULL}
                                   
                                   return(tmp_dist)
                                 })

    # print(colnames(env_dist_df))
    site_comp_df <- site_comp_df |>
      dplyr::left_join(env_dist_df |>
        dplyr::rename(Env_dist = Distance) |>
        dplyr::select(-metric_nm),
                        by = c("SSBS1", "SSBS2", "SSBS_contrast", "SS"))
  }
  


  # %%%%%%
  # TODO:
  # Custom distance calcs...
  if ("custom_dist" %in% metrics){ 

    ss_un <- unique(site_comp_df$SS)

    site_comp_df <- site_comp_df |>
      dplyr::mutate(SSBS_contrast = paste(SSBS1, SSBS2, sep = "_"))

    # Check n metrics, calc, return... 
    for (i in 1:length(custom_metric_specs)){
      tmp_spec <- custom_metric_specs[[i]]
      tmp_nm <- names(custom_metric_specs)[i]

      tmp_custom_df <- purrr::map_dfr(.x = ss_un,
                                 .f = function(ss = .x,
                                               inp_data = site_df[,c("SS", "SSBS", tmp_spec$vars)] |> 
                                                 stats::na.omit(),
                                               metr_nm = tmp_spec$dist_metr, 
                                               contrasts = site_comp_df$SSBS_contrast){
                                   
                                   tmp_df <- inp_data |> 
                                     dplyr::filter(SS == ss) |> 
                                     dplyr::select(-SS) |> 
                                     tibble::column_to_rownames("SSBS")
                                   
                                   # print(tmp_env)
                                   # Need to check that this would work even if
                                   # no env data for a whole ss.. 
                                   if (nrow(tmp_df) > 1){
                                     tmp_dist <- hib_calc_veg(tmp_df,
                                                              metr_nm,
                                                              F,
                                                              F,
                                                              F)
                                     
                                     tmp_dist <- tmp_dist |>
                                       dplyr::rename(SSBS_contrast = contrast) |>
                                       dplyr::filter(SSBS_contrast %in% contrasts)
                                     tmp_dist$SS <- ss
                                   } else {tmp_dist <- NULL}
                                   
                                   return(tmp_dist)
                                 })

    
      site_comp_df <- site_comp_df |>
        dplyr::left_join(tmp_custom_df |>
          dplyr::rename_with(~c(tmp_nm), c(Distance)) |>
          # dplyr::rename(Env_dist = Distance) |>
          dplyr::select(-metric_nm),
                          by = c("SSBS1", "SSBS2", "SSBS_contrast", "SS"))

    }
  }

  # By SS?, or over all sites... - use hib_calc_veg
  # veg_df <- purrr::map_dfr(.x = studies,
  #                                .f = function(ss = .x, 
  #                                              spp_data = spp_tmp,
  #                                              comp_df = site_comp_df,
  #                                              metric_df = veg_metrics_df){
                                   
  #                                  # print( spp_data |>
  #                                  #   # filter out the SS that matches the SS we're interested in
  #                                  #   dplyr::filter(SS == ss) |>
  #                                  #   # pull out the site name, species name and abundance information
  #                                  #   dplyr::select(SSBS, Taxon_name_entered, Measurement) |>
  #                                  #   # pivot so that each column is a species and each row is a site
  #                                  #   tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement,
  #                                  #                      values_fill = 0) )
  #                                    # TODO: fix values_fill bug...
  #                                    # tidyr::pivot_wider(names_from = Taxon_name_entered, 
  #                                    #                    values_from = Measurement) )


  #                                  spp_site_data <- spp_data |>
  #                                    # filter out the SS that matches the SS we're interested in
  #                                    dplyr::filter(SS == ss) |>
  #                                    # pull out the site name, species name and abundance information
  #                                    dplyr::select(SSBS, Taxon_name_entered, Measurement) |>
  #                                    # pivot so that each column is a species and each row is a site
  #                                    tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement,
  #                                                       values_fill = 0) |>
                                     
  #                                    # set the rownames to the SSBS and then remove that column
  #                                    tibble::column_to_rownames("SSBS")
                                   
  #                                  # TODO: check at least 2 species are present (i.e 2 cols or more)
  #                                  # TODO: return max_dist value - if NA, set rescale to F...
  #                                  # Depending on metric_vec values, calc metrics of interest
  #                                  # %%%
  #                                  out_df <- purrr::pmap(.l = metric_df,
  #                                                           .f = function(metric,
  #                                                                         binary,
  #                                                                         rel_abund,
  #                                                                         rescale,
  #                                                                         spp_data = spp_site_data){
  #                                                             # print(metric)
  #                                                             hib_calc_veg(spp_data, 
  #                                                                          metric,
  #                                                                          binary,
  #                                                                          rel_abund,
  #                                                                          rescale)
  #                                                           }) |>
  #                                    dplyr::bind_rows()
  #                                  # %%%
                                   
  #                                  out_df <- out_df |>
  #                                    dplyr::filter(contrast %in% comp_df$contrast)
                                   
  #                                  out_df <- tidyr::pivot_wider(out_df,
  #                                                               names_from = metric_nm, 
  #                                                               values_from = Distance)
  #                                  out_df$SS <- ss
                                   
  #                                  return(out_df)
  #                                })
    
  #   veg_df <- veg_df |>
  #     dplyr::select(-contrast)

  # %%%%%%


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
