#' Calculate total abundance
#' 
#' @description 
#' This function calculates total abundance for each site measuring this metric
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_ta(obj)
#' } 
#' @export
hib_calc_ta <- function(obj){
  message(crayon::green("Calculating total abundance...\n"))
  # TODO:
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(obj$species, "SS"),
    assertthat::has_name(obj$species, "SSBS"),
    assertthat::has_name(obj$species, "Measurement") || assertthat::has_name(obj$species, "Effort_corrected_measurement"),
    
    assertthat::has_name(obj$species, "Diversity_metric_type")
  )
  
  # Check if effort correction applied...
  if (obj$processing_flags$measures_effort_corrected == T){
    abund_df <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance") %>%
      dplyr::group_by(SSBS, SS) %>%
      dplyr::summarise(Total_abundance = sum(Effort_corrected_measurement)) %>%
      dplyr::ungroup()
  }
  if (obj$processing_flags$measures_effort_corrected == F){
    abund_df <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance") %>%
      dplyr::group_by(SSBS, SS) %>%
      dplyr::summarise(Total_abundance = sum(Measurement),
                       N_spp_SBSS = dplyr::n_distinct(Best_guess_binomial)) %>%
      dplyr::ungroup()
  }
  
  # Re-scale abundance to 0-1 per study (SS)
  abund_df <- abund_df %>% 
    dplyr::group_by(SS) %>%
    dplyr::mutate(Rescaled_total_abundance = Total_abundance / max(Total_abundance)) %>%
    dplyr::ungroup()
  # bind to obj$sites
  # obj$sites <- dplyr::left_join(obj$sites,
  #                               abund_df,
  #                               by = "SSBS")
  
  # Get N species used
  n_spp <- obj$species %>%
    dplyr::filter(Diversity_metric_type == "Abundance") %>%
    dplyr::select(Best_guess_binomial) %>%
    dplyr::n_distinct()
  
  # Store in obj$diversity
  out_tib <- dplyr::tibble(Metric = "Total_abundance",
                           Data = list(abund_df),
                           N_species_total = n_spp,
                           Flag = paste0("Effort corrected: ", obj$processing_flags$measures_effort_corrected))
  
  div_empty <- is.null(obj$diversity)
  
  if (!div_empty){
    obj$diversity <- dplyr::bind_rows(obj$diversity, out_tib)
  }  
  if (div_empty){
    obj$diversity <- out_tib
  }

  return(obj)
}



#' Calculate species richness
#' 
#' @description 
#' This function calculates species richness for each site
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_sr(obj)
#' } 
#' @export
hib_calc_sr <- function(obj){
  message(crayon::green("Calculating species richness...\n"))
  # TODO:
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(obj$species, "SSBS"),
    assertthat::has_name(obj$species, "Measurement") || assertthat::has_name(obj$species, "Effort_corrected_measurement"),
    
    assertthat::has_name(obj$species, "Diversity_metric_type")
  )
  
  # Check if effort correction applied...
  if (obj$processing_flags$measures_effort_corrected == T){
    sr_df1 <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance" |
                      Diversity_metric_type == "Occurence") %>%
      dplyr::group_by(SSBS) %>%
      dplyr::summarise(Species_richness = sum(Effort_corrected_measurement > 0))
    # And direct measures of spp richness
    sr_df2 <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Species richness") %>%
      dplyr::group_by(SSBS) %>%
      dplyr::summarise(Species_richness = sum(Effort_corrected_measurement))
    
    sr_df <- dplyr::bind_rows(sr_df1, sr_df2)
  }
  if (obj$processing_flags$measures_effort_corrected == F){
    sr_df1 <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance" | 
                      Diversity_metric_type == "Occurence") %>%
      dplyr::group_by(SSBS) %>%
      dplyr::summarise(Species_richness = sum(Measurement > 0))
    # And direct measures of spp richness
    sr_df2 <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Species richness") %>%
      dplyr::group_by(SSBS) %>%
      dplyr::summarise(Species_richness = sum(Measurement))
    
    sr_df <- dplyr::bind_rows(sr_df1, sr_df2)
  }
  
  # bind to obj$sites
  # obj$sites <- dplyr::left_join(obj$sites,
  #                               sr_df,
  #                               by = "SSBS")
  
  # Get N species used
  n_spp <- obj$species %>%
    dplyr::filter(Diversity_metric_type == "Abundance"  | 
                    Diversity_metric_type == "Occurence"  | 
                    Diversity_metric_type == "Species richness") %>%
    dplyr::select(Best_guess_binomial) %>%
    dplyr::n_distinct()
  
  
  # Store in obj$diversity
  out_tib <- dplyr::tibble(Metric = "Species_richness",
                           Data = list(sr_df),
                           N_species_total = n_spp,
                           Flag = paste0("Effort corrected: ", obj$processing_flags$measures_effort_corrected))
  
  div_empty <- is.null(obj$diversity)
  
  if (!div_empty){
    obj$diversity <- dplyr::bind_rows(obj$diversity, out_tib)
  }  
  if (div_empty){
    obj$diversity <- out_tib
  }
  
  return(obj)
}


#' Calculate simpson's diversity
#' 
#' @description 
#' This function calculates simpson's diversity per site
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_simpson(obj)
#' } 
#' @export
hib_calc_simpson <- function(obj){
  message(crayon::green("Calculating Simpson's diversity ...\n"))
  # TODO:
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(obj$species, "SSBS"),
    assertthat::has_name(obj$species, "Measurement") || assertthat::has_name(obj$species, "Effort_corrected_measurement"),
    
    assertthat::has_name(obj$species, "Diversity_metric_type")
  )
  
  # Check if effort correction applied...
  if (obj$processing_flags$measures_effort_corrected == T){
    #.... 
    simp_df <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance") %>%
      dplyr::group_by(SSBS) %>%
      dplyr::summarise(Simpson_diversity = 1 - sum( ( Effort_corrected_measurement / sum(Effort_corrected_measurement) )^2 ))
  }
  if (obj$processing_flags$measures_effort_corrected == F){
    #....
    simp_df <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance") %>%
      dplyr::group_by(SSBS) %>%
      dplyr::summarise(Simpson_diversity = 1 - sum( ( Measurement / sum(Measurement) )^2 ))
  }
  
  # bind to obj$sites
  # obj$sites <- dplyr::left_join(obj$sites,
  #                               simp_df,
  #                               by = "SSBS")
  
  # Get N species used
  n_spp <- obj$species %>%
    dplyr::filter(Diversity_metric_type == "Abundance") %>%
    dplyr::select(Best_guess_binomial) %>%
    dplyr::n_distinct()
  
  # Store in obj$diversity
  out_tib <- dplyr::tibble(Metric = "Simpson_diversity",
                           Data = list(simp_df),
                           N_species_total = n_spp,
                           Flag = paste0("Effort corrected: ", obj$processing_flags$measures_effort_corrected))
  
  div_empty <- is.null(obj$diversity)
  
  if (!div_empty){
    obj$diversity <- dplyr::bind_rows(obj$diversity, out_tib)
  }  
  if (div_empty){
    obj$diversity <- out_tib
  }
  
  return(obj)
}

# TODO:
# Chao-1 richness
# Rarefied richenss
# %%%%%
# From https://github.com/timnewbold/predicts-demo/blob/master/predictsFunctions/R/SiteMetrics.R
# 
# .SiteChao <- function(diversity, site.abundance) {
#   
#   cat("Computing Chao\n")
#   # Compute Chao only for studies where Diversity_metric_is_suitable_for_Chao 
#   # is TRUE and all measurements within the study are integers.
#   
#   # http://stackoverflow.com/questions/3476782/how-to-check-if-the-number-is-integer
#   # TODO all.equal
#   study.chao.suitable <- 
#     tapply(diversity$Measurement[diversity$Diversity_metric_is_suitable_for_Chao],
#            diversity$SS[diversity$Diversity_metric_is_suitable_for_Chao], 
#            function(m) all(floor(m)==m))
#   
#   # Studies for which Diversity_metric_is_suitable_for_Chao is FALSE will have 
#   # a value of NA in study.chao.suitable
#   study.chao.suitable[is.na(study.chao.suitable)] <- FALSE
#   site.chao.suitable <- study.chao.suitable[tapply(diversity$SS, diversity$SSS, unique)]
#   
#   diversity$Is_Chao_suitable <- site.chao.suitable[diversity$SSS]
#   
#   # Only abundance metrics should be marked as suitable for Chao 
#   stopifnot(!any(site.chao.suitable & !site.abundance))
#   chao <- rep(NA, length(site.abundance))
#   chao[site.chao.suitable] <- 
#     tapply(diversity$Measurement[diversity$Is_Chao_suitable], 
#            droplevels(diversity$SSS[diversity$Is_Chao_suitable]), 
#            function(m) sum(m>0) + (((sum(m==1) * (sum(m==1)-1)) / (2*(sum(m==2)+1)))))
#   return (chao)
# }

#' Calculate Chao-1 species richness
#' 
#' @description 
#' This function calculates Chao-1 species richness per site
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_chao1(obj)
#' } 
#' @export
hib_calc_chao1 <- function(obj){
  message(crayon::green("Calculating Chao1 species richness ...\n"))
  # TODO:
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(obj$species, "SSBS"),
    assertthat::has_name(obj$species, "SS"),
    assertthat::has_name(obj$species, "Measurement") || assertthat::has_name(obj$species, "Effort_corrected_measurement"),
    
    assertthat::has_name(obj$species, "Diversity_metric_type")
  )
  
  # Check if effort correction applied...
  if (obj$processing_flags$measures_effort_corrected == T){
    #.... 
    # Filter to suitable sites...
    # Abundance, Chao-suitable (all integers...)
    chao1_df <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance" &
                      Diversity_metric_is_suitable_for_Chao) %>%
      dplyr::group_by(SS) %>%
      dplyr::filter(all(floor(Effort_corrected_measurement) == Effort_corrected_measurement)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(SSBS) %>%
      dplyr::summarise(Chao1 = sum(Effort_corrected_measurement > 0) + (((sum(Effort_corrected_measurement == 1) * (sum(Effort_corrected_measurement == 1) - 1)) / (2 * (sum(Effort_corrected_measurement == 2) + 1)))) ) 
    # More efficient to use eqn...
    # unname(vegan::estimateR(Effort_corrected_measurement)["S.chao1"]) 
    # 
  }
  if (obj$processing_flags$measures_effort_corrected == F){
    chao1_df <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance" &
                      Diversity_metric_is_suitable_for_Chao) %>%
      dplyr::group_by(SS) %>%
      dplyr::filter(all(floor(Measurement) == Measurement)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(SSBS) %>%
      dplyr::summarise(Chao1 = sum(Measurement > 0) + (((sum(Measurement == 1) * (sum(Measurement == 1) - 1)) / (2 * (sum(Measurement == 2) + 1)))) ) 
    # unname(vegan::estimateR(Measurement)["S.chao1"]) )
  }
  
  # bind to obj$sites
  # obj$sites <- dplyr::left_join(obj$sites,
  #                               simp_df,
  #                               by = "SSBS")
  n_spp <- obj$species %>%
    dplyr::filter(Diversity_metric_type == "Abundance" &
                    Diversity_metric_is_suitable_for_Chao) %>%
    dplyr::group_by(SS) %>%
    dplyr::filter(all(floor(Measurement) == Measurement)) %>%
    dplyr::ungroup() %>%
    dplyr::select(Best_guess_binomial) %>%
    dplyr::n_distinct()
  
  # Store in obj$diversity
  out_tib <- dplyr::tibble(Metric = "Chao1_diversity",
                           Data = list(chao1_df),
                           N_species_total = n_spp,
                           Flag = paste0("Effort corrected: ", obj$processing_flags$measures_effort_corrected))
  
  div_empty <- is.null(obj$diversity)
  
  if (!div_empty){
    obj$diversity <- dplyr::bind_rows(obj$diversity, out_tib)
  }  
  if (div_empty){
    obj$diversity <- out_tib
  }
  
  return(obj)
}

# Rarefied richness

# .SiteRarefiedRichness<-function(diversity, site.abundance){
#   
#   cat("Computing Rarefied Species Richness\n")
#   # Compute rarefied richness only for studies where Diversity_metric_is_suitable_for_Chao 
#   # is TRUE and all measurements within the study are integers.
#   
#   # http://stackoverflow.com/questions/3476782/how-to-check-if-the-number-is-integer
#   # TODO all.equal
#   study.rsr.suitable <- 
#     tapply(diversity$Measurement[diversity$Diversity_metric_is_suitable_for_Chao],
#            diversity$SS[diversity$Diversity_metric_is_suitable_for_Chao], 
#            function(m) all(floor(m)==m))
#   
#   # Studies for which Diversity_metric_is_suitable_for_Chao is FALSE will have 
#   # a value of NA in study.rsr.suitable
#   study.rsr.suitable[is.na(study.rsr.suitable)] <- FALSE
#   site.rsr.suitable <- study.rsr.suitable[tapply(diversity$SS, diversity$SSS, unique)]
#   
#   diversity$Is_RSR_suitable <- site.rsr.suitable[diversity$SSS]
#   
#   # Only abundance metrics should be marked as suitable for rarefied richness 
#   stopifnot(!any(site.rsr.suitable & !site.abundance))
#   
#   rsrich <- rep(NA, length(site.abundance))
#   
#   if (any(diversity$Is_RSR_suitable)){
#     siteTotalAbund<-aggregate(Measurement~SSS,data=droplevels(
#       diversity[diversity$Is_RSR_suitable,]),FUN=sum)
#     siteTotalAbund$SS<-diversity$SS[match(siteTotalAbund$SSS,diversity$SSS)]
#     studyMinSiteAbund<-aggregate(Measurement~SS,data=siteTotalAbund,FUN=
#                                    function(x) min(x[x>0]))
#     diversity$studyMinSiteAbund<-studyMinSiteAbund$Measurement[match(
#       diversity$SS,studyMinSiteAbund$SS)]
#     
#     sample.n <- split(diversity$studyMinSiteAbund[diversity$Is_RSR_suitable], 
#                       droplevels(diversity$SSS[diversity$Is_RSR_suitable]))
#     values <- split(diversity$Measurement[diversity$Is_RSR_suitable],
#                     droplevels(diversity$SSS[diversity$Is_RSR_suitable]))
#     
#     rsr<-matrix(nrow=length(which(site.rsr.suitable)),ncol=1000)
#     for (i in 1:1000){
#       rsr[,i]<-mapply(function(vals,samp){
#         sp<-as.character(1:length(vals))
#         ind<-rep(sp,vals)
#         if (length(ind)==0){
#           sr<-0
#         } else {
#           rare.sp<-sample(ind,samp[1])
#           sr<-length(table(rare.sp))
#         }
#         return(sr)
#       },values,sample.n)
#     }
#     
#     rsrich[site.rsr.suitable]<-apply(rsr,1,mean)
#     
#   }
# 
#   return(rsrich)
#   
# }
# %%%%%

#' Calculate rarefied species richness
#' 
#' @description 
#' This function calculates rarefied species richness per site
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_rsr(obj)
#' } 
#' @export
hib_calc_rsr <- function(obj){
  message(crayon::green("Calculating rarefied species richness ...\n"))
  # TODO:
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(obj$species, "SSBS"),
    assertthat::has_name(obj$species, "SS"),
    assertthat::has_name(obj$species, "Measurement") || assertthat::has_name(obj$species, "Effort_corrected_measurement"),
    
    assertthat::has_name(obj$species, "Diversity_metric_type")
  )
  
  rarefy_comm <- function(abund_distr, N){
    comm <- rep(letters[1:length(abund_distr)], abund_distr)
    samp <- sample(comm, N)
    sr <- length(table(samp))
    return(sr)
  }
  
  # Check if effort correction applied...
  if (obj$processing_flags$measures_effort_corrected == T){
    #.... 
    # Filter to suitable sites...
    # Abundance, Chao-suitable (all integers...)
    rsr_df <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance" &
                      Diversity_metric_is_suitable_for_Chao) %>%
      dplyr::group_by(SSBS) %>%
      dplyr::mutate(TA = sum(Effort_corrected_measurement)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(SS) %>%
      dplyr::mutate(min_abund = min(TA)) %>%
      dplyr::filter(all(floor(Effort_corrected_measurement) == Effort_corrected_measurement)) %>%
      dplyr::ungroup() 
    
    # Need to repeatedly sample from each SSBS abundance distr - 1000x
    # expand spp abund to long vect
    # sample based on min_abund
    # calc sr per SSBS/repeat
    rsr_df <- purrr::map_dfr(seq(1,1000),
                             function(x, df){
                               df %>%
                                 dplyr::group_by(SSBS) %>%
                                 dplyr::summarise(sr = rarefy_comm(Effort_corrected_measurement, 
                                                                   min_abund[1]),
                                                  min_abund = min_abund[1],
                                                  iter = x)
                             },
                             df = rsr_df) %>%
      dplyr::group_by(SSBS) %>%
      dplyr::summarise(mean_rsr = mean(sr),
                       min_abund = min_abund[1])
  }
  
  if (obj$processing_flags$measures_effort_corrected == F){
    # Filter to suitable sites...
    # Abundance, Chao-suitable (all integers...)
    rsr_df <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance" &
                      Diversity_metric_is_suitable_for_Chao) %>%
      dplyr::group_by(SSBS) %>%
      dplyr::mutate(TA = sum(Measurement)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(SS) %>%
      dplyr::mutate(min_abund = min(TA)) %>%
      dplyr::filter(all(floor(Measurement) == Measurement)) %>%
      dplyr::ungroup() 
    
    # Need to repeatedly sample from each SSBS abundance distr - 1000x
    # expand spp abund to long vect
    # sample based on min_abund
    # calc sr per SSBS/repeat
    rsr_df <- purrr::map_dfr(seq(1,1000),
                             function(x, df){
                               df %>%
                                 dplyr::group_by(SSBS) %>%
                                 dplyr::summarise(sr = rarefy_comm(Measurement, 
                                                                   min_abund[1]),
                                                  min_abund = min_abund[1],
                                                  iter = x)
                             },
                             df = rsr_df) %>%
      dplyr::group_by(SSBS) %>%
      dplyr::summarise(mean_rsr = mean(sr),
                       min_abund = min_abund[1])
  }
  
  # bind to obj$sites
  # obj$sites <- dplyr::left_join(obj$sites,
  #                               simp_df,
  #                               by = "SSBS")
  
  n_spp <- obj$species %>%
    dplyr::filter(Diversity_metric_type == "Abundance" &
                    Diversity_metric_is_suitable_for_Chao) %>%
    dplyr::group_by(SS) %>%
    dplyr::filter(all(floor(Measurement) == Measurement)) %>%
    dplyr::ungroup() %>%
    dplyr::select(Best_guess_binomial) %>%
    dplyr::n_distinct()
  # Store in obj$diversity
  out_tib <- dplyr::tibble(Metric = "Rarefied_species_richness",
                           Data = list(rsr_df),
                           N_species_total = n_spp,
                           Flag = paste0("Effort corrected: ", obj$processing_flags$measures_effort_corrected))
  
  div_empty <- is.null(obj$diversity)
  
  if (!div_empty){
    obj$diversity <- dplyr::bind_rows(obj$diversity, out_tib)
  }  
  if (div_empty){
    obj$diversity <- out_tib
  }
  
  return(obj)

}


#' Calculate similarity metrics
#' 
#' @description 
#' This function calculates either Sorensen, Jaccard, Bray-curtis or Ruzika 
#' similarity metrics between sites, using the betapart package
#' @details 
#' Expects a dataframe/tibble with species abundance records at two sites 
#' - sites as rows, species as columns
#' @param spp_site_data A dataframe/tibble object.
#' @param metric A [`character`] of either sorensen, jaccard, bray, ruzicka, hypervolume
#' @param use_relative_abundance A [`logical`] value - whether relative abundance 
#' values should be used; applicable for "bray", "ruzicka", and "hypervolume"
#' @param rescale_hypervolume  A [`logical`] value - whether hypervolume should be
#' divided by 0.5 (theoretical maximum volume for 2 communities)
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   metric_val <- hib_calc_sim_(spp_site_data, "sorensen")
#' } 
#' @export
hib_calc_sim_ <- function(spp_site_data, metric, 
                          use_relative_abundance,
                          rescale_hypervolume){
  
  assertthat::assert_that(
    is.character(metric)
    # assertthat::has_name(spp_site_data, "SSBS"),
    # assertthat::has_name(spp_site_data, "Taxon_name_entered"),
    # assertthat::has_name(spp_site_data, "Measurement")
  )
  
  metric <- match.arg(metric, c("jaccard", "sorensen", "bray", "ruzicka",
                                "hypervolume", "chao-jaccard", "chao-sorensen"), 
                      several.ok = FALSE)
  
  # if (use_relative_abund == T){
  #   spp_site_data <- spp_site_data/rowSums(spp_site_data)
  # }
  # print(metric)
  # print(spp_site_data)
  if (metric %in% c("jaccard", "sorensen")){
    spp_site_data[spp_site_data > 0] <- 1
    # print(spp_site_data)
    out <- 1 - unlist(betapart::beta.pair(spp_site_data, 
                                          index.family = metric))
  } else if (metric %in% c("bray", "ruzicka")){
    if (use_relative_abundance){
      spp_site_data <- spp_site_data/rowSums(spp_site_data)
    }
    out <- 1 - unlist(betapart::beta.pair.abund(spp_site_data, 
                                                index.family = metric))
  }
  # TODO
  # add ChaoJaccard, ChaoSorensen
  else if (metric %in% c("chao-jaccard", "chao-sorensen")){
    # transpose data for use in SpadeR
    spp_site_data <- as.data.frame(t(spp_site_data))
    # print(as.data.frame(spp_site_data))
    # print(class(as.data.frame(spp_site_data)))
    out <- SpadeR::SimilarityPair(X = spp_site_data, datatype = "abundance")
    # Get the metric of interest...and rename columns
    if (metric == "chao-jaccard"){
      out <- as.data.frame.list(out$estimated_relative["ChaoJaccard",])
      colnames(out) <- c("ChaoJaccard", "ChaoJaccard_se", 
                         "ChaoJaccard_lo95", "ChaoJaccard_hi95")
    }
    if (metric == "chao-sorensen"){
      out <- as.data.frame.list(out$estimated_relative["ChaoSorensen",])
      colnames(out) <- c("ChaoSorensen", "ChaoSorensen_se", 
                         "ChaoSorensen_lo95", "ChaoSorensen_hi95")
    }
    # Output is in similarity form (1 == comms are equal)
    return(out)
  }
  
  # Add in hypervolume approach
  else if (metric == "hypervolume"){
    
    # transpose site_spp_data for use in convhulln (M rows of N dimensions)
    spp_site_data <- t(spp_site_data)
    
    if (use_relative_abundance){
      spp_site_data <- spp_site_data/colSums(spp_site_data)
    }
    # add origin to spp_site_data
    spp_site_data <- rbind(spp_site_data,
                           rep(0, ncol(spp_site_data)))
    # vol <- geometry::convhulln(spp_site_data,
    #                            output.options = T)$vol
    
    # %%%
    vol <- tryCatch(
      {
        geometry::convhulln(spp_site_data, 
                            output.options = TRUE)$vol
      },
      error = function(e) {
        0
      }
    )
    # %%%
    if (rescale_hypervolume == T){
      vol <- vol/0.5
    }
    out <- data.frame("Hypervolume" = vol)
    return (out)
  }

  out <- as.data.frame.list(out)
  # print(out)
  # if one of the sites doesn't have any individuals/species in it
  # i.e. the row sum is 0
  if (sum(rowSums(spp_site_data) == 0, na.rm = TRUE) == 1){
    # then the similarity between sites should be 0
    out[1,] <- 0
    # if both sites have no individuals/species
  } else if (sum(rowSums(spp_site_data) == 0, na.rm = TRUE) == 2){
    # then class the similarity as NA
    out[1,] <- NA
  }
  
  # Rename columns
  if (metric == "sorensen"){
    out <- out %>%
      dplyr::rename(Simpson = beta.sim,
                    Simpson_nested = beta.sne,
                    Sorensen = beta.sor)
  }
  if (metric == "jaccard"){
    out <- out %>%
      dplyr::rename(Jaccard_turnover = beta.jtu,
                    Jaccard_nested = beta.jne,
                    Jaccard = beta.jac)
  }
  if (metric == "bray"){
    out <- out %>%
      dplyr::rename(Balanced_bray_curtis = beta.bray.bal,
                    Gradient_bray_curtis = beta.bray.gra,
                    Bray_curtis = beta.bray)
  }
  if (metric == "ruzicka"){
    out <- out %>%
      dplyr::rename(Balanced_ruzicka = beta.ruz.bal,
                    Gradient_ruzicka = beta.ruz.gra,
                    Ruzicka = beta.ruz)
  }
  
  return (out)
}

#' Calculate compositional similarity (wrapper)
#' 
#' @description 
#' This function calculates compositional similarity between sites within a study
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @param sim_metric A similarity metric, any/all of c("balanced bray-curtis", "bray-curtis", ...)
#' @param primary_minimal_baseline A [`boolean`] value, determines whether 
#' compositional similarity is calculated for all site combinations per study, or
#' only relative to a primary minimal baseline
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_sim(obj, c("balanced bray-curtis"))
#' } 
#' @export
hib_calc_sim <- function(obj, 
                         sim_metric = "bray", 
                         primary_minimal_baseline = T,
                         use_relative_abundance = F,
                         rescale_hypervolume = T,
                         rescale_euclidean_distance = T){
  message(crayon::green("Calculating compositional similarity...\n"))
  # TODO:
  # Assertion checks
  assertthat::assert_that(
    is.character(sim_metric),
    assertthat::has_name(obj$species, "SSBS"),
    assertthat::has_name(obj$species, "SS"),
    assertthat::has_name(obj$species, "Measurement"),
    
    assertthat::has_name(obj$species, "Predominant_land_use"),
    assertthat::has_name(obj$species, "Taxon_name_entered"),
    
    assertthat::has_name(obj$sites, "SSBS"),
    assertthat::has_name(obj$sites, "SS")
    # ?Effort_correct_measurement?
    # Use_intensity
  )
  
  sim_metric <- match.arg(sim_metric, c("jaccard", "sorensen", "bray", "ruzicka",
                                        "euclidean", "hypervolume",
                                        "chao-jaccard", "chao-sorensen"), several.ok = TRUE)
  
  # Workflow from : https://adrianadepalma.github.io/BII_tutorial/bii_example.html
  # Only use known land-use sites
  # Only use abundance metric data
  # Only use studies (SS) with single sampling effort - use Measurement...
  # Only use studies with >1 species
  # Only use studies with primary minimal land-use
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
  
  if (primary_minimal_baseline == T){
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
    
  } else if (primary_minimal_baseline == F){
    site_tmp <- obj$sites %>%
      dplyr::filter(SSBS %in% unique(spp_tmp$SSBS)) %>%
      # left_join(unique(spp_tmp[, c("SSBS", "prim_min")])) %>%
      droplevels()
  }
  
  # TODO: don't what this to affect chaoJaccard/chaoSorensen
  # if (use_relative_abundance){
  #   spp_tmp <- spp_tmp %>% 
  #     group_by(SS, SSBS) %>% 
  #     mutate(Measurement = Measurement/sum(Measurement)) %>%
  #     dplyr::ungroup()
  #   }
  
  
  n_spp <- spp_tmp %>%
    dplyr::select(Best_guess_binomial) %>%
    dplyr::n_distinct()
    
  # Need to get combinations of data... - s1, s2, ...
  # %%%
  # From https://adrianadepalma.github.io/BII_tutorial/bii_example.html
  
  # get a vector of each study to loop over
  studies <- unique(as.character(site_tmp$SS))
  
  # print(studies)
  # %>%
  #   dplyr::distinct(SS) %>%
  #   dplyr::pull()
  
  site_comp_df <- purrr::map_dfr(.x = studies, .f = function(x){
    
    # filter out the given study
     #%>%
      # simplify the data so we only have one row for each site
      # dplyr::distinct(SSBS, .keep_all = TRUE)
    
    # print(site_data)
    
    # pull out the sites that are Primary minimal (we only want to use comparisons with the baseline)
    if (primary_minimal_baseline == T){
      site_data <- dplyr::filter(site_tmp, SS == x) %>%
        # pull out the SSBS and LandUse information
        dplyr::select(SSBS, prim_min)
      
      baseline_sites <- site_data %>%
        dplyr::filter(prim_min == 1) %>%
        dplyr::pull(SSBS)
    } else if (primary_minimal_baseline == F){
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
      dplyr::rename(s1 = Var1, s2 = Var2) %>%
      
      # remove the comparisons where the same site is being compared to itself
      dplyr::filter(s1 != s2) %>%
      
      # make the values characters rather than factors
      dplyr::mutate(s1 = as.character(s1),
                    s2 = as.character(s2),
                    
                    # add the full name
                    contrast = paste(s1, s2, sep = "_"),
                    
                    # add the study id
                    SS = as.character(x))

    return(site_comparisons)
  })
  # print(head(site_comp_df))
  message(crayon::green(paste0(nrow(site_comp_df), 
                               " site comparisons to make...\n"))
          )
  # # Need furrr
  # future::plan(future::multisession, workers = 4)
  # 
  # # I'm using map2 (because there are two arguments I'm passing through - s1 and s2)
  # # and the map2_dbl because the output I want is a vector of numbers (double rather than integer format)
  # # so this function is going to go through each s1 and s2 in turn
  # # and pass them into the get_bray function
  # 
  # # # %%%%%
  # # # RC: this bit is relevant to any/all metrics
  # # sp_data <- data %>%
  # #   
  # #   # filter out the SSBS that matches the pair of sites we're interested in
  # #   dplyr::filter(SSBS %in% c(s1, s2)) %>%
  # #   
  # #   # pull out the site name, species name and abundance information
  # #   dplyr::select(SSBS, Taxon_name_entered, Measurement) %>%
  # # # %%%%# 
  #   
  # sim_df <- furrr::future_map2_dfr(.x = site_comp_df$s1,
  #                         .y = site_comp_df$s2,
  #                         .f = function(s1 = .x, s2 = .y, 
  #                                       spp_data = spp_tmp,
  #                                       metric_vec = sim_metric){
  #                           require(hibiiscus)
  #                           require(betapart)
  #                           spp_site_data <- spp_data %>%
  #                             # filter out the SSBS that matches the pair of sites we're interested in
  #                             dplyr::filter(SSBS %in% c(s1, s2)) %>%
  #                             # pull out the site name, species name and abundance information
  #                             dplyr::select(SSBS, Taxon_name_entered, Measurement)
  #                           
  #                           out_df <- data.frame(s1 = s1, 
  #                                                s2 = s2)
  #                           
  #                           # Depending on sim_metric values, calc metrics of comp-sim
  #                           if (sum(grepl("bray", metric_vec)) > 0 ){
  #                             bray_metr_vec <- metric_vec[grepl("bray", metric_vec)]
  #                             tmp_bray <- hib_calc_sim_bray(spp_site_data, bray_metr_vec)
  #                             out_df <- dplyr::bind_cols(out_df, tmp_bray)
  #                           }
  #                           # Other sim metrics - sorensen?...
  #                           
  #                           
  #                           return(out_df)
  #                         }
  #                         )
  # 
  # # stop running things in parallel for now
  # future::plan(future::sequential)
  # # %%%
  
  
  # %%%
  # Non-parallel option - V slow for even a 1000 sites...
  sim_df <- NULL
  
  non_euc_metr <- sim_metric[sim_metric %in% c("jaccard", "sorensen", "bray", "ruzicka",
                                               "hypervolume", "chao-jaccard", "chao-sorensen")]
  # TODO
  # add ChaoJaccard, ChaoSorensen
  
  if (length(non_euc_metr)>0){
    
    sim_df <- purrr::map2_dfr(.x = site_comp_df$s1,
                              .y = site_comp_df$s2,
                              .f = function(s1 = .x, s2 = .y, 
                                            spp_data = spp_tmp,
                                            metric_vec = non_euc_metr){
                                # require(hibiiscus)
                                # require(betapart)
                                spp_site_data <- spp_data %>%
                                  # filter out the SSBS that matches the pair of sites we're interested in
                                  dplyr::filter(SSBS %in% c(s1, s2)) %>%
                                  # pull out the site name, species name and abundance information
                                  dplyr::select(SSBS, Taxon_name_entered, Measurement) %>%
                                  # pivot so that each column is a species and each row is a site
                                  tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement) %>%
                                  # set the rownames to the SSBS and then remove that column
                                  tibble::column_to_rownames("SSBS")
                                
                                out_df <- data.frame(s1 = s1, 
                                                     s2 = s2)
                                
                                # print(paste0(s1, s2, sep = "_"))
                                
                                # Depending on sim_metric values, calc metrics of comp-sim
                                out_met <- purrr::map_dfc(.x = metric_vec,
                                                          .f = function(metric = .x,
                                                                        spp_data){
                                                            hib_calc_sim_(spp_data, metric,
                                                                          use_relative_abundance,
                                                                          rescale_hypervolume)
                                                          },
                                                          spp_data = spp_site_data)
                                out_df <- out_df %>% bind_cols(out_met)
                                
                                return(out_df)
                              }
    )
  }
  
  # TODO:
  # Note: for euclidean distance, can just take subset of all site combs per study...
  # Use relative abundance, rescale - add as flags for euclidean distance...
  # But need to know which are prim_min... and only keep those... makes sense to just incl in above
  # But v inefficient...
  # Better to calc all distance within study, add  site idx to long version
  euc_sim_df <- NULL
  if ("euclidean" %in% sim_metric){
    
    if (use_relative_abundance){
      spp_tmp_euc <- spp_tmp %>%
        group_by(SS, SSBS) %>%
        mutate(Measurement = Measurement/sum(Measurement)) %>%
        dplyr::ungroup()
      }

    euc_sim_df <- purrr::map_dfr(.x = studies,
                                 .f = function(ss = .x, 
                                               spp_data = spp_tmp_euc,
                                               comp_df = site_comp_df){
                                   spp_site_data <- spp_data %>%
                                     # filter out the SS that matches the SS we're interested in
                                     dplyr::filter(SS == ss) %>%
                                     # pull out the site name, species name and abundance information
                                     dplyr::select(SSBS, Taxon_name_entered, Measurement) %>%
                                     # pivot so that each column is a species and each row is a site
                                     tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement) %>%
                                     # set the rownames to the SSBS and then remove that column
                                     tibble::column_to_rownames("SSBS")
                                   
                                   # if (use_relative_abund == T){
                                   #   spp_site_data <- spp_site_data/rowSums(spp_site_data)
                                   # }
                                   
                                   # print(head(spp_site_data))
                                   dist_matr <- as.matrix(dist(spp_site_data))
                                   
                                   rownames(dist_matr) <- rownames(spp_site_data)
                                   colnames(dist_matr) <- rownames(spp_site_data)
                                   
                                   # print(head(dist_matr))
                                   
                                   dist_df <- as.data.frame(dist_matr)
                                   dist_df$s1 <- rownames(dist_df)
                                   dist_df <- tidyr::pivot_longer(dist_df, cols = !s1, 
                                                                  names_to = "s2", 
                                                                  values_to = "Euclidean_distance")
                                   
                                   dist_df <- dist_df %>%
                                     dplyr::filter(s1 != s2) %>% 
                                     dplyr::mutate(contrast = paste(s1, s2, sep = "_")) %>%
                                     dplyr::filter(contrast %in% comp_df$contrast)
                                   # print(head(dist_df))
                                   # dist_df$SS <- ss
                                   # Subset to comps of interest
                                   # dist_df <- dplyr::filter(dist_df,
                                   #                          contrast %in% comp_df$contrast)
                                   
                                   if (rescale_euclidean_distance == T){
                                     dist_df$Euclidean_distance <- dist_df$Euclidean_distance / sqrt(2)
                                   }
                                   
                                   return(dist_df)
                                 })
  }
  
  # Combine metrics into single df
  if ((!is.null(sim_df)) & (!is.null(euc_sim_df))){
    sim_df <- left_join(sim_df, euc_sim_df[,c("s1", "s2", "Euclidean_distance")])
  }
  if (is.null(sim_df) & (!is.null(euc_sim_df))){
    sim_df <- euc_sim_df[,c("s1", "s2", "Euclidean_distance")]
  }
  
  
  # %%%
  # calc metrics of interest...
  # can be v slow....
  # need to multithread...
  
  # TODO: other options to add...
  
  
  # Need to make output longer...
  # row in obj$diversity per metric
  metrics_out <- colnames(sim_df)
  metrics_out <- metrics_out[!metrics_out %in% c("s1", "s2", 
                                                 "ChaoJaccard_se", 
                                                 "ChaoJaccard_lo95", "ChaoJaccard_hi95",
                                                 "ChaoSorensen_se", 
                                                 "ChaoSorensen_lo95", "ChaoSorensen_hi95")]
  
  metrics_out_ls <- purrr::map(metrics_out,
                               function(x,
                                        df = sim_df){
                                 if (x %in% c("ChaoJaccard", "ChaoSorensen")){
                                   return(df[,c("s1", "s2", x,
                                                paste0(x, "_se"),
                                                paste0(x, "_lo95"),
                                                paste0(x, "_hi95"))])
                                 }
                                 return(df[,c("s1", "s2", x)])
                               })
  
  out_tib <- dplyr::tibble(Metric = metrics_out,
                           Data = metrics_out_ls,
                           N_species_total = n_spp,
                           Flag = paste0("Primary minimal baseline: ", primary_minimal_baseline,
                                         "; Using relative abundance: ", use_relative_abundance,
                                         "; Rescaling hypervolumes: ", rescale_hypervolume,
                                         "; Rescaling euclidean distance: ", rescale_euclidean_distance))
  
  # diversity$Metric[x] as <metric>
  # diversity$Data[[x]] should be [s1, s2, <metric>]
  # ....
  
  div_empty <- is.null(obj$diversity)
  
  if (!div_empty){
    obj$diversity <- dplyr::bind_rows(obj$diversity, out_tib)
  }  
  if (div_empty){
    obj$diversity <- out_tib
  }
  
  return(obj)
}


