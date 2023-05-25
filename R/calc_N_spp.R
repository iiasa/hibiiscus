#' Calculate number of species contributing to a model - internal function
#' 
#' @description 
#' This function calculates the number of species contributing to a model fit based
#' on the metric (model response) of interest
#' @details 
#' Expects a [`Database`] object.
#' @param obj A [`Database`] object.
#' @param mod_df A [`data.frame`] object created as part of data preparation for
#' model fitting.
#' @param mod_resp A [`character`] object indicating the modelled response variable.
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   n_spp <- calc_N_spp(obj, mod_df, mod_resp)
#' } 
#' @noRd
calc_N_spp <- function(obj,
                       mod_df,
                       mod_resp){
  
  # First, for within site metrics
  if (mod_resp %in% c("Total_abundance", "Rescaled_total_abundance", 
                      "Simpson_diversity")){
    n_spp <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance") %>%
      dplyr::filter(SSBS %in% mod_df$SSBS) %>%
      dplyr::select(Best_guess_binomial) %>%
      dplyr::n_distinct()
  }
  
  if (mod_resp == "Species_richness"){
    n_spp <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance"  | 
                      Diversity_metric_type == "Occurence"  | 
                      Diversity_metric_type == "Species richness") %>%
      dplyr::filter(SSBS %in% mod_df$SSBS) %>%
      dplyr::select(Best_guess_binomial) %>%
      dplyr::n_distinct()
  }
  if (mod_resp %in% c("Chao1_diversity", "Rarefied_species_richness")){
    n_spp <- obj$species %>%
      dplyr::filter(Diversity_metric_type == "Abundance" &
                      Diversity_metric_is_suitable_for_Chao) %>%
      dplyr::group_by(SS) %>%
      dplyr::filter(all(floor(Measurement) == Measurement)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(SSBS %in% mod_df$SSBS) %>%
      dplyr::select(Best_guess_binomial) %>%
      dplyr::n_distinct()
  }
  
  # Or for between site metrics
  if (mod_resp %in% tools::toTitleCase(c("manhattan", "euclidean", "canberra", "clark", "bray", 
                                         "kulczynski", "jaccard", "gower", "altGower", 
                                         "morisita", "horn", "mountford", "raup", "binomial", 
                                         "chao", "cao", "mahalanobis", "chisq", "chord", 
                                         "hellinger", "aitchison", "robust.aitchison",
                                         
                                         "hypervolumne", "cosine",
                                         "baselga_bray_balanced", "baselga_bray_gradient",
                                         "baselga_ruzicka_balanced", "baselga_ruzicka_gradient"))){
    
    n_spp <- obj$species %>%
      dplyr::filter(!is.na(Predominant_land_use) &
                      Diversity_metric_type == "Abundance") %>%
      dplyr::group_by(SS) %>%
      dplyr::mutate(n_sample_effort = n_distinct(Sampling_effort),
                    # Taxon_name_entered v best_guess_binomial ??
                    n_species = dplyr::n_distinct(Taxon_name_entered)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n_sample_effort == 1 &
                      n_species > 1) %>%
      dplyr::filter(SSBS %in% unique(c(mod_df$SSBS1, mod_df$SSBS2))) %>%
      dplyr::select(Best_guess_binomial) %>%
      dplyr::n_distinct()
  }
  
  return(n_spp)
  
}

