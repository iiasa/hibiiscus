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
      dplyr::summarise(Total_abundance = sum(Measurement)#,
                       # N_spp_SBSS = dplyr::n_distinct(Best_guess_binomial)
                       ) %>%
      dplyr::ungroup()
  }
  
  # Re-scale abundance to 0-1 per study (SS)
  abund_df <- abund_df %>% 
    dplyr::group_by(SS) %>%
    dplyr::mutate(Rescaled_total_abundance = Total_abundance / max(Total_abundance)) %>%
    dplyr::ungroup()

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
  
  return(out_tib)
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
  
  return(out_tib)
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
  message(crayon::green("Calculating within-site Simpson's diversity ...\n"))
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

  return(out_tib)
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
      dplyr::summarise(Chao1_diversity = sum(Effort_corrected_measurement > 0) + (((sum(Effort_corrected_measurement == 1) * (sum(Effort_corrected_measurement == 1) - 1)) / (2 * (sum(Effort_corrected_measurement == 2) + 1)))) ) 
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
      dplyr::summarise(Chao1_diversity = sum(Measurement > 0) + (((sum(Measurement == 1) * (sum(Measurement == 1) - 1)) / (2 * (sum(Measurement == 2) + 1)))) ) 
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
  
  # TODO:
  return(out_tib)
  
  # div_empty <- is.null(obj$diversity)
  # 
  # if (!div_empty){
  #   obj$diversity <- dplyr::bind_rows(obj$diversity, out_tib)
  # }  
  # if (div_empty){
  #   obj$diversity <- out_tib
  # }
  # 
  # return(obj)
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
      dplyr::summarise(Rarefied_species_richness = mean(sr),
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
      dplyr::summarise(Rarefied_species_richness = mean(sr),
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
  
  return(out_tib)

}


#' Calculate within-site diversity metrics
#' 
#' @description 
#' This function calculates diversity per site, wrapper for the above functions
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @param metrics A [`character`] vector, any of "Total_abundance", "Species_richness", 
#' "Chao1_diversity",  "Simpson_diversity", "Rarefied_species_richness"
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_div_within(obj, c("Total_abundance", "Species_richness"))
#' } 
#' @export
hib_calc_div_within <- function(obj, metrics){
  # Assertions
  assertthat::assert_that(
    assertthat::has_name(obj$species, "SSBS"),
    assertthat::has_name(obj$species, "SS"),
    assertthat::has_name(obj$species, "Measurement") || assertthat::has_name(obj$species, "Effort_corrected_measurement"),
    assertthat::has_name(obj$species, "Diversity_metric_type"),
    is.character(metrics)
  )
  
  # Match supplied metrics to options.
  metrics <- match.arg(metrics,
                       c("Total_abundance", "Species_richness", "Chao1_diversity",  
                         "Simpson_diversity", "Rarefied_species_richness"),
                       several.ok = T)
  
  message(crayon::green("Calculating within-site diversity metrics ...\n"))
  
  out_tib <- purrr::map(.x = metrics, 
                       .f = function(metric = .x, db_obj = obj){
                         if (metric == "Total_abundance"){
                           return(hib_calc_ta(db_obj))
                         }
                         if (metric == "Species_richness"){
                           return(hib_calc_sr(db_obj))
                         }
                         if (metric == "Chao1_diversity"){
                           return(hib_calc_chao1(db_obj))
                         }
                         if (metric == "Simpson_diversity"){
                           return(hib_calc_simpson(db_obj))
                         }
                         if (metric == "Rarefied_species_richness"){
                           return(hib_calc_rsr(db_obj))
                         }
                         }
                       ) %>% bind_rows()
  
  div_empty <- is.null(obj$diversity_within)
  
  if (!div_empty){
    obj$diversity_within <- dplyr::bind_rows(obj$diversity_within, out_tib)
  }  
  if (div_empty){
    obj$diversity_within <- out_tib
  }
  
  return(obj)
}
# %%%%%


# TODO: similarity v within-site metrics should be stored separately..
# diversity_within - focussed around SSBS
# diversity_between - focussed around SSBS_1 v SSBS_2 - key cols - SSBS_contrast, SS, <Metric>



# TODO: vegan R package can calc many "distance" metrics in an approach analagous
# to the dist() function, may be much quicker, includes Horn-Morista...
# Takes a site (rows) by spp (columns) matrix as input

# vegan metrics:
# c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", 
# "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", 
# "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", "robust.aitchison")

# conflicts between current metrics and vegan metrics...
# c("jaccard", "bray")

# what about sorensen - what is this?? - binary bray...
# ruzika - equivalent to quantitative jaccard??

# also need to record if one site has 0 total abundance...

# use_relative_abundance should only apply to euclidean and hypervolume.. - maybe bray..
# all other metrics specifically include some sort of scaling


# TODO:
# %%%%%%%%%%%%%

#' Calculate vegan-based dissimilarity metrics
#' 
#' @description 
#' This function calculates vegan-based dissimilarity metrics between sites
#' @details 
#' Expects a dataframe/tibble with species abundance records at two sites 
#' - sites as rows, species as columns
#' @param spp_site_data A dataframe/tibble object.
#' @param metric A [`character`] - see vegdist
#' @param binary A [`logical`] value - should distance be based on presence/absence
#' @param rel_abund A [`logical`] value - whether relative abundance values should 
#' be used; applicable for euclidean ...
#' @param rescale A [`logical`] value - should distance be rescaled relative to
#' theoretical max (if applicable and not 1)
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   metric_val <- hib_calc_veg(spp_site_data, "euclidean")
#' } 
#' @export
hib_calc_veg <- function(spp_site_data, 
                         metric,
                         binary,
                         rel_abund,
                         rescale
                         ){
  
  assertthat::assert_that(
    is.character(metric),
    is.logical(binary) || is.null(binary),
    is.logical(rel_abund) || is.null(rel_abund),
    is.logical(rescale) || is.null(rescale)
    # assertthat::has_name(spp_site_data, "SSBS"),
    # assertthat::has_name(spp_site_data, "Taxon_name_entered"),
    # assertthat::has_name(spp_site_data, "Measurement")
  )
  # print(use_relative_abundance)
  
  metric <- match.arg(metric, 
                      c("manhattan", "euclidean", "canberra", "clark", "bray", 
                        "kulczynski", "jaccard", "gower", "altGower", 
                        "morisita", "horn", "mountford", "raup", "binomial", 
                        "chao", "cao", "mahalanobis", "chisq", "chord", 
                        "hellinger", "aitchison", "robust.aitchison"), 
                      several.ok = FALSE)
  
  
  # if ((rel_abund) & (metric == "euclidean")){
  #   row_abund <- rowSums(spp_site_data)
  #   row_abund[row_abund==0] <- 1
  #   spp_site_data <- spp_site_data/row_abund
  # }
  if (rel_abund){
    row_abund <- rowSums(spp_site_data)
    row_abund[row_abund==0] <- 1
    spp_site_data <- spp_site_data/row_abund
  }
  
  
  dist_obj <- vegan::vegdist(spp_site_data,
                             method = metric,
                             binary = binary)
  
  # Extract distance matrix
  dist_matr <- as.matrix(dist_obj)
  max_dist <- attr(dist_obj, "maxdist")
  rm(dist_obj)
  
  rownames(dist_matr) <- rownames(spp_site_data)
  colnames(dist_matr) <- rownames(spp_site_data)
  
  # Convert to long
  dist_df <- as.data.frame(dist_matr)
  dist_df$SSBS1 <- rownames(dist_df)
  dist_df <- tidyr::pivot_longer(dist_df, cols = !SSBS1, 
                                 names_to = "SSBS2", 
                                 values_to = "Distance")
  
  dist_df <- dist_df %>%
    dplyr::filter(SSBS1 != SSBS2) %>% 
    dplyr::mutate(metric_nm = tools::toTitleCase(metric),
                  contrast = paste(SSBS1, SSBS2, sep = "_"))
  
  if (rescale & !is.na(max_dist)){
    dist_df$Distance <- dist_df$Distance / max_dist
  }
  # Ideally, if rescale is T and max_dist is NA, edit rescale to be F in output
  
  
  # TODO: if as_sim, convert to similarity score (not possible for unscaled distances...)
  
  return(dist_df)
}



#' Calculate non-vegan-based dissimilarity metrics
#' 
#' @description 
#' This function calculates non-vegan-based dissimilarity metrics between sites
#' @details 
#' Expects a dataframe/tibble with species abundance records at two sites 
#' - sites as rows, species as columns
#' @param spp_site_data A dataframe/tibble object.
#' @param metric A [`character`] - cosine, hypervolume
#' @param binary A [`logical`] value - should distance be based on presence/absence
#' @param rel_abund A [`logical`] value - whether relative abundance values should 
#' be used; applicable for ...hypervolume
#' @param rescale A [`logical`] value - should distance be rescaled relative to
#' theoretical max (if applicable and not 1)
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   metric_val <- hib_calc_nonveg(spp_site_data, "cosine")
#' } 
#' @export
hib_calc_nonveg <- function(spp_site_data, 
                            metric,
                            binary,
                            rel_abund,
                            rescale
                            ){
  
  assertthat::assert_that(
    is.character(metric),
    is.logical(binary) || is.null(binary),
    is.logical(rel_abund) || is.null(rel_abund),
    is.logical(rescale) || is.null(rescale)
    # assertthat::has_name(spp_site_data, "SSBS"),
    # assertthat::has_name(spp_site_data, "Taxon_name_entered"),
    # assertthat::has_name(spp_site_data, "Measurement")
  )
  # print(use_relative_abundance)
  metric <- match.arg(metric, 
                      c("hypervolumne", "cosine",
                        "baselga_bray_balanced", "baselga_bray_gradient", 
                        "baselga_ruzicka_balanced", "baselga_ruzicka_gradient"),
                      # "baselga_bray","baselga_ruzicka",
                      # "baselga_jaccard", "baselga_sorensen") 
                      several.ok = FALSE)
  
  # Add in hypervolume approach
  if (metric == "hypervolume"){
    
    if (rel_abund){
      row_abund <- rowSums(spp_site_data)
      row_abund[row_abund == 0] <- 1
      spp_site_data <- spp_site_data/row_abund
      # spp_site_data <- spp_site_data/colSums(spp_site_data)
    }
    if (binary){
      spp_site_data[spp_site_data > 0] <- 1
    }
    
    # transpose site_spp_data for use in convhulln (M rows of N dimensions)
    spp_site_data <- t(spp_site_data)
    
    # add origin to spp_site_data
    spp_site_data <- rbind(spp_site_data,
                           rep(0, ncol(spp_site_data)))

    vol <- tryCatch(
      {
        geometry::convhulln(spp_site_data, 
                            output.options = TRUE)$vol
      },
      error = function(e) {
        0
      }
    )
    
    if (rescale & rel_abund){
      vol <- vol/0.5
    }
    
    out <- data.frame("Hypervolume" = vol)
    # return (out)
  }
  
  else if (metric == "cosine"){
    x <- as.matrix(spp_site_data)
    # N rows == 2
    x1 <- x[1,]
    x2 <- x[2,]
    cosine <- sum( x1 * x2 ) / ( sqrt( sum( x1^2 )) * sqrt( sum(x2^2) ) )
    
    # if one of the sites doesn't have any individuals/species in it
    # i.e. the row sum is 0
    if (sum(rowSums(spp_site_data) == 0, na.rm = TRUE) == 1){
      # then the similarity between sites should be 0
      cosine <- 0
      # if both sites have no individuals/species
    } else if (sum(rowSums(spp_site_data) == 0, na.rm = TRUE) == 2){
      # then class the similarity as NA
      cosine <- NA
    }
    out <- data.frame("Cosine" = cosine)
    # return (out)
  }
  
  
  # else if (metric %in% c("jaccard", "sorensen")){
  #   # "baselga_jaccard", "baselga_sorensen"
  #   spp_site_data[spp_site_data > 0] <- 1
  #   # print(spp_site_data)
  #   out <- 1 - unlist(betapart::beta.pair(spp_site_data, 
  #                                         index.family = metric))
  # } 
  else if (metric %in% c("baselga_bray_balanced", "baselga_bray_gradient",
                         "baselga_ruzicka_balanced", "baselga_ruzicka_gradient")){
    # "baselga_bray","baselga_ruzicka",
    
    # (metric %in% c("baselga_bray_balanced", "baselga_bray_gradient",
                  # "baselga_ruzicka_balanced", "baselga_ruzicka_gradient"))
    
    # stop("In development")
    
    if (rel_abund){
      # rowSums could be 0 (zero)!!! - adjust equal 1 if it does...
      row_abund <- rowSums(spp_site_data)
      row_abund[row_abund == 0] <- 1
      spp_site_data <- spp_site_data/row_abund
    }
    if (binary){
      spp_site_data[spp_site_data > 0] <- 1
    }
    
    if (metric %in% c("baselga_bray_balanced", "baselga_bray_gradient")){
      out <- unlist(betapart::beta.pair.abund(spp_site_data, 
                                              index.family = "bray"))
    }
    if (metric %in% c("baselga_ruzicka_balanced", "baselga_ruzicka_gradient")){
      out <- unlist(betapart::beta.pair.abund(spp_site_data, 
                                              index.family = "ruzicka"))
    }
    
    out <- as.data.frame.list(out)
    # print(out)
    # if one of the sites doesn't have any individuals/species in it
    # i.e. the row sum is 0
    # TODO: but is the below true for the gradient/turnover components calculated
    # by betapart...
    # if (sum(rowSums(spp_site_data) == 0, na.rm = TRUE) == 1){
    #   # then the similarity between sites should be 0
    #   out[1,] <- 0
    #   # if both sites have no individuals/species
    # } else if (sum(rowSums(spp_site_data) == 0, na.rm = TRUE) == 2){
    #   # then class the similarity as NA
    #   out[1,] <- NA
    # }
    
    # Rename columns
    if (metric %in% c("baselga_bray_balanced", "baselga_bray_gradient")){
      out <- out %>%
        dplyr::rename(baselga_bray_balanced = beta.bray.bal,
                      baselga_bray_gradient = beta.bray.gra)
    }
    if (metric %in% c("baselga_ruzicka_balanced", "baselga_ruzicka_gradient")){
      out <- out %>%
        dplyr::rename(baselga_ruzicka_balanced = beta.ruz.bal,
                      baselga_ruzicka_gradient = beta.ruz.gra)
    }
    
    out <- out %>%
      dplyr::select(!!metric)
    
    colnames(out) <- tools::toTitleCase(colnames(out))
  }
  
  # TODO:
  # fix above...
  
  # TODO: if as_sim, convert to similarity score (not possible for unscaled distances...)
  
  
  return(out)
}


#' Calculate compositional dissimilarity (wrapper) between sites
#' 
#' @description 
#' This function calculates compositional dissimilarity between sites within a study
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @param metrics A [`character`] vector of dissimilarity metrics, any/all of ...
#' @param binary A [`logical`] vector applicable for vegan-based metrics, see [`vegdist`], 
#' for non-vegan metrics, this will be ignored
#' @param rel_abund A [`logical`] vector determining whether relative abundance 
#' should be used to calculate the dissimilarity metric, only applicable for some metrics, 
#' e.g., euclidean, hypervolume - Note: many metrics include rescaling to 'relative abundances'
#' as part of the calculation
#' @param rescale A [`logical`] vector determining if dissimilarity measure should 
#' be rescaled relative to its maximum value - applies to metrics that are not 
#' on a 0-1 scale by default - hypervolume, euclidean (if relative abundance), chord, hellinger
#' @param baseline_lui A [`character`] or NULL, indicating the baseline LUI against
#' which compositional dis/similarity is calculated per study
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_div_between(obj, 
#'                               c("bray", "jaccard"),
#'                               c(F, T),
#'                               c(T, F),
#'                               c(F, F),
#'                               "Primary vegetation_Minimal use")
#' } 
#' @export
hib_calc_div_between <- function(obj, 
                                 metrics,
                                 binary = F,
                                 rel_abund = T,
                                 rescale = T,
                                 # primin_baseline = T
                                 baseline_lui = "Primary vegetation_Minimal use"
                                 ){

  # Assertions
  assertthat::assert_that(
    
    inherits(obj, "Database"),
    # input checks
    is.character(metrics),
    is.logical(binary),
    is.logical(rel_abund),
    is.logical(rescale),
    is.character(baseline_lui) || is.null(baseline_lui),
    
    assertthat::has_name(obj$species, "SSBS"),
    assertthat::has_name(obj$species, "SS"),
    assertthat::has_name(obj$species, "Measurement"),
    assertthat::has_name(obj$species, "Diversity_metric_type"),
    
    assertthat::has_name(obj$species, "Predominant_land_use"),
    assertthat::has_name(obj$species, "Use_intensity"),
    assertthat::has_name(obj$species, "LUI"),
    assertthat::has_name(obj$species, "Taxon_name_entered"),
    
    assertthat::has_name(obj$sites, "SSBS"),
    assertthat::has_name(obj$sites, "SS")
  )
  
  message(crayon::green("Preparing site-contrasts to calculate compositional dis/similarity...\n"))
  
  veg_metrics <- c("manhattan", "euclidean", "canberra", "clark", "bray", 
                   "kulczynski", "jaccard", "gower", "altGower", 
                   "morisita", "horn", "mountford", "raup", "binomial", 
                   "chao", "cao", "mahalanobis", "chisq", "chord", 
                   "hellinger", "aitchison", "robust.aitchison")
                   
                   # Non-Vegan  
  non_veg_metrics <- c("hypervolumne", "cosine",
                       "baselga_bray_balanced", "baselga_bray_gradient",
                       "baselga_ruzicka_balanced", "baselga_ruzicka_gradient")
  # "baselga_bray","baselga_ruzicka",
  # "baselga_jaccard", "baselga_sorensen")
  
  metrics_ <- c(veg_metrics, non_veg_metrics)
  
  # first, combine metrics, binary, rel_abund and rescale into a df
  metrics_spec_df <- data.frame(metric    = metrics,
                                binary    = binary,
                                rel_abund = rel_abund,
                                rescale   = rescale)
  
  matchable_metrics_idx <- !is.na(pmatch(metrics_spec_df$metric,
                                         metrics))
  
  # subset to matchable metrics
  metrics_spec_df <- metrics_spec_df[matchable_metrics_idx,]
  
  # run match.arg on matchable metrics
  metrics_spec_df$metric <- match.arg(metrics_spec_df$metric, 
                                       metrics_,
                                       several.ok = T)
  
  # split to veg and non_veg spec dfs
  veg_metrics_df <- subset(metrics_spec_df,
                           metric %in% veg_metrics)
  non_veg_metrics_df <- subset(metrics_spec_df,
                               metric %in% non_veg_metrics)
  
  # Get site comparisons
  spp_tmp <- obj$species
  
  spp_tmp <- spp_tmp %>%
    dplyr::filter(!is.na(Predominant_land_use) &
                    Diversity_metric_type == "Abundance") %>%
    dplyr::group_by(SS) %>%
    dplyr::mutate(n_sample_effort = n_distinct(Sampling_effort),
           # Taxon_name_entered v best_guess_binomial ??
           n_species = dplyr::n_distinct(Taxon_name_entered)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n_sample_effort == 1 &
                    n_species > 1)
  
  # if (primin_baseline == T){
  #   # TODO: should this be ammended to primary forest...
  #   spp_tmp$prim_min <- as.numeric((grepl("primary", spp_tmp$Predominant_land_use, ignore.case = T) &
  #                                     grepl("minimal", spp_tmp$Use_intensity, ignore.case = T)))
  #   spp_tmp <- spp_tmp %>%
  #     dplyr::group_by(SS) %>%
  #     dplyr::mutate(n_primin_records = sum(prim_min)) %>%
  #     dplyr::filter(n_primin_records > 0) %>%
  #     dplyr::ungroup()
  #   
  #   site_tmp <- obj$sites %>%
  #     dplyr::filter(SSBS %in% unique(spp_tmp$SSBS)) %>%
  #     left_join(unique(spp_tmp[, c("SSBS", "prim_min")])) %>%
  #     droplevels()
  #   
  # } else if (primin_baseline == F){
  #   site_tmp <- obj$sites %>%
  #     dplyr::filter(SSBS %in% unique(spp_tmp$SSBS)) %>%
  #     droplevels()
  # }
  if (!is.null(baseline_lui)){
    spp_tmp$base_lui <- spp_tmp$LUI == baseline_lui
    
    spp_tmp <- spp_tmp %>%
      dplyr::group_by(SS) %>%
      dplyr::mutate(n_base_lui = sum(base_lui)) %>%
      dplyr::filter(n_base_lui > 0) %>%
      dplyr::ungroup()
    
    site_tmp <- obj$sites %>%
      dplyr::filter(SSBS %in% unique(spp_tmp$SSBS)) %>%
      dplyr::left_join(unique(spp_tmp[, c("SSBS", "base_lui")])) %>%
      droplevels()
    
  } else if (is.null(baseline_lui)){
    site_tmp <- obj$sites %>%
      dplyr::filter(SSBS %in% unique(spp_tmp$SSBS)) %>%
      droplevels()
  }
  
  n_spp <- spp_tmp %>%
    dplyr::select(Best_guess_binomial) %>%
    dplyr::n_distinct()
  
  # Need to get combinations of data... - s1, s2, ...
  # %%%
  # From https://adrianadepalma.github.io/BII_tutorial/bii_example.html
  
  # get a vector of each study to loop over
  studies <- unique(as.character(site_tmp$SS))
  
  site_comp_df <- purrr::map_dfr(.x = studies, .f = function(x){
    
    # pull out the sites that are the baseline lui (we only want to use comparisons with the baseline)
    if (!is.null(baseline_lui)){
      site_data <- dplyr::filter(site_tmp, SS == x) %>%
        # pull out the SSBS and LandUse information
        dplyr::select(SSBS, base_lui)
      
      baseline_sites <- site_data %>%
        dplyr::filter(base_lui == 1) %>%
        dplyr::pull(SSBS)
    } else if (is.null(baseline_lui)){
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
  
  message(crayon::green(paste0(nrow(site_comp_df), 
                               " site comparisons to make...\n"))
          )
  
  message(crayon::green("Calculating compositional dis/similarity...\n"))
  
  # do vegan-based metrics...
  # %%%
  veg_df <- NULL
  
  if (nrow(veg_metrics_df) > 0){
    message(crayon::green("Calculating vegan-based dissimilarity metrics...\n"))
    
    veg_df <- purrr::map_dfr(.x = studies,
                                 .f = function(ss = .x, 
                                               spp_data = spp_tmp,
                                               comp_df = site_comp_df,
                                               metric_df = veg_metrics_df){
                                   # TODO: add other processing flags too!!
                                   # binary, rel_abund, rescale
                                   
                                   spp_site_data <- spp_data %>%
                                     # filter out the SS that matches the SS we're interested in
                                     dplyr::filter(SS == ss) %>%
                                     # pull out the site name, species name and abundance information
                                     dplyr::select(SSBS, Taxon_name_entered, Measurement) %>%
                                     # pivot so that each column is a species and each row is a site
                                     tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement,
                                                        values_fill = 0) %>%
                                     # set the rownames to the SSBS and then remove that column
                                     tibble::column_to_rownames("SSBS")
                                   
                                   # Depending on metric_vec values, calc metrics of interest
                                   # %%%
                                   out_df <- purrr::pmap(.l = metric_df,
                                                            .f = function(metric,
                                                                          binary,
                                                                          rel_abund,
                                                                          rescale,
                                                                          spp_data = spp_site_data){
                                                              # print(metric)
                                                              hib_calc_veg(spp_data, 
                                                                           metric,
                                                                           binary,
                                                                           rel_abund,
                                                                           rescale)
                                                            }) %>%
                                     dplyr::bind_rows()
                                   # %%%
                                   
                                   out_df <- out_df %>%
                                     dplyr::filter(contrast %in% comp_df$contrast)
                                   
                                   out_df <- tidyr::pivot_wider(out_df,
                                                                names_from = metric_nm, 
                                                                values_from = Distance)
                                   out_df$SS <- ss
                                   
                                   return(out_df)
                                 })
    
    veg_df <- veg_df %>%
      dplyr::select(-contrast)
  }
  # %%%
  
  # do other metrics...
  non_veg_df <- NULL
  
  if (nrow(non_veg_metrics_df) > 0){
    message(crayon::green("Calculating non-vegan dissimilarity metrics...\n"))
    
    # stop("TODO")
    non_veg_df <- purrr::map2_dfr(.x = site_comp_df$s1,
                              .y = site_comp_df$s2,
                              .f = function(s1 = .x, 
                                            s2 = .y, 
                                            spp_data = spp_tmp,
                                            metric_d = non_veg_metrics_df){

                                spp_site_data <- spp_data %>%
                                  # filter out the SSBS that matches the pair of sites we're interested in
                                  dplyr::filter(SSBS %in% c(s1, s2)) %>%
                                  # pull out the site name, species name and abundance information
                                  dplyr::select(SSBS, Taxon_name_entered, Measurement) %>%
                                  # pivot so that each column is a species and each row is a site
                                  tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement,
                                                     values_fill = 0) %>%
                                  # set the rownames to the SSBS and then remove that column
                                  tibble::column_to_rownames("SSBS")
                                
                                ss <- spp_data %>%
                                  # filter out the SSBS that matches the pair of sites we're interested in
                                  dplyr::filter(SSBS %in% c(s1, s2)) %>%
                                  dplyr::pull(SS) %>%
                                  unique()
                                
                                out_df <- data.frame(SSBS1 = s1, 
                                                     SSBS2 = s2,
                                                     SS    = ss)
                                
                                # Depending on sim_metric values, calc metrics of comp-sim
                                # non_veg_metrics_df
                                out_met <- purrr::pmap(.l = non_veg_metrics_df,
                                                          .f = function(metric,
                                                                        binary,
                                                                        rel_abund,
                                                                        rescale,
                                                                        spp_data = spp_site_data){
                                                            hib_calc_nonveg(spp_data, 
                                                                            metric,
                                                                            binary,
                                                                            rel_abund,
                                                                            rescale)
                                                          }) %>%
                                  dplyr::bind_cols()
                                
                                out_df <- out_df %>% bind_cols(out_met)
                                return(out_df)
                              }
    )
  }
  
  # combine veg and non_veg
  # s1, s2, ...
  if ((!is.null(veg_df)) & (!is.null(non_veg_df))){
    veg_df <- left_join(veg_df, non_veg_df) 
  }
  if (is.null(veg_df) & (!is.null(non_veg_df))){
    veg_df <- non_veg_df
  }
  
  # convert to long tibble...
  # add processing flags per metric - 
  
  # Calc is_empty flags once and add to each sub df?
  empty_sites1 <- spp_tmp %>%
    dplyr::group_by(SSBS) %>%
    dplyr::select(SSBS, Taxon_name_entered, Measurement) %>%
    # pivot so that each column is a species and each row is a site
    tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement,
                       values_fill = 0) %>%
    tibble::column_to_rownames("SSBS") %>%
    dplyr::mutate(SSBS1_is_empty = (rowSums(.) == 0)) %>%
    tibble::rownames_to_column("SSBS1") %>%
    dplyr::select(SSBS1, SSBS1_is_empty)
  
  empty_sites2 <- empty_sites1 %>%
    dplyr::rename(SSBS2 = SSBS1,
                  SSBS2_is_empty = SSBS1_is_empty)
  
    
  veg_df <- veg_df %>%
    dplyr::left_join(empty_sites1) %>%
    dplyr::left_join(empty_sites2)
  # 
  metrics_out <- colnames(veg_df)
  metrics_out <- metrics_out[!metrics_out %in% c("SSBS1", "SSBS2", "SS",
                                                 "SSBS1_is_empty", "SSBS2_is_empty")] 

  metrics_out_ls <- purrr::map(metrics_out,
                               function(x,
                                        df = veg_df){
                                 return(df[,c("SSBS1", "SSBS2", "SS", x,
                                              "SSBS1_is_empty", "SSBS2_is_empty")])
                               })
  
  out_tib <- dplyr::tibble(Metric = metrics_out,
                           Data = metrics_out_ls,
                           N_species_total = n_spp
                           )
  # Add Flag based on metric and specific flags...
  # need to match sentence-case metric names to lower-case metric names...
  metrics_spec_df <- metrics_spec_df %>%
    dplyr::mutate(Metric = tools::toTitleCase(metric))
  
  metrics_spec_df <- metrics_spec_df %>%
    dplyr::mutate(baseline_lui = baseline_lui)
  if (is.null(baseline_lui)){
    metrics_spec_df <- metrics_spec_df %>%
      dplyr::mutate(baseline_lui = NA)
  }
  
  # Improve column names...
  metrics_spec_df_ <- metrics_spec_df[,c("Metric", "binary", "rel_abund", "rescale",
                                         "baseline_lui")]
  colnames(metrics_spec_df_) <- c("Metric", "Binary", "Relative_abundance", "Rescale",
                                  "Baseline_LUI")
  out_tib <- out_tib %>%
    dplyr::left_join(metrics_spec_df_)
  
  
  div_empty <- is.null(obj$diversity_between)
  
  if (!div_empty){
    obj$diversity_between <- dplyr::bind_rows(obj$diversity_between, out_tib)
  }  
  if (div_empty){
    obj$diversity_between <- out_tib
  }
  
  return(obj)
}

