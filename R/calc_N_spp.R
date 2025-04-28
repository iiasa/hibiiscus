#' calc_N_spp
#' 
#' @description 
#' Calculate number of species contributing to a model (internal function).
#' 
#' @details 
#' This function calculates the number of species contributing to a model fit 
#' based on the metric (model response) of interest. Intended for internal use.
#' Calls \code{get_n_spp_taxa}.
#' 
#' @param obj A HibDatabase object.
#' @param mod_df A data.frame object created as part of data preparation for
#' model fitting.
#' @param mod_resp A character object indicating the modelled response variable.
#' 
#' @return list
#' 
#' @keywords database, model, processing
#' 
#' @seealso
#' [get_n_spp_taxa()]
#' 
#' @examples 
#' \dontrun{
#'   n_spp <- calc_N_spp(obj, mod_df, mod_resp)
#' } 
#' 
#' @noRd
calc_N_spp <- function(obj,
                       mod_df,
                       mod_resp){
  # TODO: add binary for between site responses...??

  assertthat::assert_that(
    inherits(obj, "HibDatabase"),
    is.data.frame(mod_df),
    is.character(mod_resp)
    )
  
  # Check if measurement or effort corrected - is this needed?
  # When would one be >0 and the other not
  if (obj$processing_flags$measures_effort_corrected == T){
    measure_var <- "Effort_corrected_measurement"
  } 
  if (obj$processing_flags$measures_effort_corrected == F){
    measure_var <- "Measurement"
  }

  if (all(c("SSBS1", "SSBS2") %in% colnames(mod_df))){
    ssbs <- unique(c(mod_df$SSBS1, mod_df$SSBS2))
  } else {
    ssbs <- unique(mod_df$SSBS)
  }

  # TODO: isn't it more efficient to just filter obj$species to SSBS/SSBS1/SSBS2 in mod_df?
  # This should be okay..
  # then calc N_spp/N_taxa
  df <- obj$species |>
    dplyr::filter(SSBS %in% ssbs)
  out <- get_n_spp_taxa(df, measure_var)

  # Surely no SSBS will have multiple types of Diversity_metric_type ?
  # Then don't even need mod_resp

  # Make more general:
  # filter df based on response...
# 
#   # First, for within site metrics
#   if (mod_resp %in% c("Total_abundance", "Rescaled_total_abundance", 
#                       "Simpson_diversity", 
#                       "Mean_abundance", "Geometric_mean_abundance",
#                       "Rescaled_mean_abundance", "Rescaled_geometric_mean_abundance",
#                       "Shannon_diversity",
#                       "HillShannon_diversity", "HillSimpson_diversity",
#                       "PIE", "sPIE", "Evenness")){
#     df <- obj$species |>
#       dplyr::filter(Diversity_metric_type == "Abundance")
#   }
#   
#   if (mod_resp == "Species_richness"){
#     df <- obj$species |>
#       dplyr::filter(Diversity_metric_type == "Abundance"  | 
#                       Diversity_metric_type == "Occurence"  | 
#                       Diversity_metric_type == "Species richness")
#   }
# 
#   if (mod_resp %in% c("Chao1_diversity", "Rarefied_species_richness")){
#     df <- obj$species |>
#       dplyr::filter(Diversity_metric_type == "Abundance" &
#                       Diversity_metric_is_suitable_for_Chao) |>
#       dplyr::group_by(SS) |>
#       dplyr::filter(all(floor(Measurement) == Measurement)) |>
#       dplyr::ungroup() 
# 
#   }
#   
#   # Or for between site metrics
#   if (mod_resp %in% tools::toTitleCase(c("manhattan", "euclidean", "canberra", "clark", "bray", 
#                                          "kulczynski", "jaccard", "gower", "altGower", 
#                                          "morisita", "horn", "mountford", "raup", "binomial", 
#                                          "chao", "cao", "mahalanobis", "chisq", "chord", 
#                                          "hellinger", "aitchison", "robust.aitchison",
#                                          
#                                          "hypervolume", "cosine", "theta",
#                                          "baselga_bray_balanced", "baselga_bray_gradient",
#                                          "baselga_ruzicka_balanced", "baselga_ruzicka_gradient",
#                                          "horn_dissim", "normalised_euclidean"
#                                          ))){
#     # TODO: check binary or not...
# 
#     df <- obj$species |>
#       dplyr::filter(!is.na(Predominant_land_use) &
#                       Diversity_metric_type == "Abundance") |>
#       dplyr::group_by(SS) |>
#       dplyr::mutate(n_sample_effort = dplyr::n_distinct(Sampling_effort),
#                     # Taxon_name_entered v best_guess_binomial ??
#                     n_species = dplyr::n_distinct(Taxon_name_entered)) |>
#       dplyr::ungroup() |>
#       dplyr::filter(n_sample_effort == 1 &
#                       n_species > 1)
#     # TODO: reconfigure - not essential to have n_sample_effort == 1
#     # %%% - Based on between-site metric calcs, note, also need binary FLAG
# 
#     # spp_tmp <- obj$species |>
#     #   dplyr::filter(#!is.na(Predominant_land_use) & # ?? - drop this
#     #                 Diversity_metric_type %in% c("Abundance",
#     #                                              "Occurrence")) |>
#       
#     #   dplyr::group_by(SS) |>
#     #   dplyr::mutate(#n_sample_effort = dplyr::n_distinct(Sampling_effort),
#     #                 # Taxon_name_entered v best_guess_binomial ??
#     #                 n_species = dplyr::n_distinct(Taxon_name_entered)) |>
#     #   dplyr::ungroup() |>
#     #   dplyr::filter(#n_sample_effort == 1 &
#     #                   n_species > 1)
#     
#     # # print(metrics_spec_df)
#     # if (all(metrics_spec_df$binary == FALSE)){
#     #   spp_tmp <- spp_tmp |>
#     #     dplyr::filter(Diversity_metric_type == "Abundance")
#     #   # print(nrow(spp_tmp))
#     #   if (nrow(spp_tmp) == 0){
#     #     stop("\nError: All metric specifications require Abundance records (binary = FALSE), but none are available.", call. = F)
#     #   }
#     # }
#     # %%%
#   }

  

  # then calc N_spp/N_taxa
  # df <- df |>
  #   dplyr::filter(SSBS %in% ssbs)
  #    
  # out <- get_n_spp_taxa(df, measure_var)

  return(out)
  
}

