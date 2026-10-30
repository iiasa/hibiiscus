# Internal functions used to estimate a range of diversity metrics


#' get_n_spp_taxa
#' 
#' @description 
#' Calculate number of species (best guess) and original taxa names.
#' 
#' @details 
#' This function calculates number of species contributing to diversity metrics.
#' It expects the relevant subset of the species data.frame from a HibDatabase.
#' Note, only species present in a site (measure > 0) contribute.
#' 
#' @param df A subset of the species data.frame from a [`HibDatabase`].
#' @param measure_var Character indicating whether to use raw ('Measurement'), 
#' or effort corrected measurements ('Effort_corrected_measurement').
#' 
#' @return list
#' 
#' @keywords database, processing
#' 
#' @examples 
#' \dontrun{
#'   out <- get_n_spp_taxa(df, "Measurement")
#' } 
#' 
#' @noRd
get_n_spp_taxa <- function(df, measure_var){
  # Assertion checks
  assertthat::assert_that(
    is.data.frame(df),
    is.character(measure_var)
  )

  assertthat::assert_that(
    assertthat::has_name(df, "Taxon_name_entered"),
    assertthat::has_name(df, "Best_guess_binomial"),
    assertthat::has_name(df, measure_var)
  )

  # Only consider present species
  df <- df |>
    dplyr::filter(!!dplyr::sym(measure_var) > 0)

  # Separately get species number based on 'Best_guess_binomial'
  n_spp <- df |>
    dplyr::filter(Best_guess_binomial != '' &
                  !is.na(Best_guess_binomial)) |>
    dplyr::select(Best_guess_binomial) |>
    dplyr::n_distinct()

  # And 'Taxon_name_entered'
  n_orig_taxa <- df |>
    dplyr::filter(Taxon_name_entered != '' &
                  !is.na(Taxon_name_entered)) |>
    dplyr::select(Taxon_name_entered) |>
    dplyr::n_distinct()

  return(list(N_spp = n_spp,
              N_taxa = n_orig_taxa))
}


#' hib_calc_ta
#' 
#' @description 
#' Calculate total abundance for each site with abundance records.
#'
#' @details 
#' Total abundance = sum(`x[i]`), where `x[i]` is the abundance of species i within
#' a given site.
#' Rescaled total abundance = `Total abundance[j,k]`/max(`Total abundance[j,k]`), where
#' the total abundance in eack site j is divided by the maximum total abundance 
#' in study k.
#' Internal function for use in `hib_calc_div_within()`.
#'
#' @param spp_df The species data.frame from a HibDatabase object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within `hib_calc_div_within()`.
#'
#' @return tibble
#' 
#' @keywords database, processing, diversity
#' 
#' @seealso
#' [hib_calc_div_within()]
#'
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_ta(df, "Measurement")
#' } 
#'
#' @noRd
hib_calc_ta <- function(spp_df, measure_var){
  
  message(crayon::white("\tTotal abundance ...\n"))
   
  # Assertion checks
  assertthat::assert_that(
    is.data.frame(spp_df),
    is.character(measure_var)
  )

  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"),
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )

  # Filter to useable records
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance") 

  abund_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::summarise(Total_abundance = sum(!!dplyr::sym(measure_var))) |>
    dplyr::ungroup() |> 
    # Re-scale abundance to 0-1 per study (SS)
    dplyr::group_by(SS) |>
    # if max Total_abundance is 0..., leave as is..
    dplyr::mutate(max_ta = max(Total_abundance),
                  Rescaled_total_abundance = dplyr::case_when(max_ta == 0 ~ Total_abundance,
                                                              T ~ Total_abundance / max_ta)) |>
    dplyr::ungroup() |>
    dplyr::select(-max_ta)
  
  # Ensure all rescaled values between 0 and 1! - avoid rounding errors
  abund_df$Rescaled_total_abundance[abund_df$Rescaled_total_abundance < 0] <- 0
  abund_df$Rescaled_total_abundance[abund_df$Rescaled_total_abundance > 1] <- 1
  
  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  # Return tibble for use in diversity_within
  out_tib <- dplyr::tibble(Metric = "Total_abundance",
                           Data = list(abund_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa)

  return(out_tib)
}



#' Calculate species richness
#' 
#' @description 
#' This function calculates species richness for each site.
#'
#' @details 
#' Species richness is the number of species (or taxonomic entities) that are 
#' present in a site.
#' Internal function for use in `hib_calc_div_within()`
#' 
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement',
#' generated within `hib_calc_div_within()`.
#' 
#' @return tibble
#' 
#' @keywords database, processing, diversity
#' 
#' @seealso
#' [hib_calc_div_within()]
#' 
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_sr(spp_df, measure_var)
#' } 
#' 
#' @noRd
hib_calc_sr <- function(spp_df, measure_var){
  message(crayon::white("\tSpecies richness ...\n"))
  
  # Assertion checks
  assertthat::assert_that(
    is.data.frame(spp_df),
    is.character(measure_var)
  )

  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), # obj$species
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance" |
                    Diversity_metric_type == "Occurence" |
                    Diversity_metric_type == "Species richness" ) 

  sr_df <- spp_df |>
    dplyr::mutate(n_sp = dplyr::case_when(Diversity_metric_type %in% c("Abundance", "Occurence") ~ as.numeric(!!dplyr::sym(measure_var) > 0),
                                          Diversity_metric_type == "Species richness" ~ !!dplyr::sym(measure_var),
                                          T ~ 0)) |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::summarise(Species_richness = sum(n_sp)) |>
    dplyr::ungroup()
  
  # Get N species used
  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  out_tib <- dplyr::tibble(Metric = "Species_richness",
                           Data = list(sr_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa)                           
  
  return(out_tib)
}


#' hib_calc_simpson
#' 
#' @description 
#' Calculate Simpson's diversity index per site.
#'
#' @details 
#' Simpson's diversity is calculated as: `1 - sum(p[i]^2)`,
#' where `p[i] = x[i] / sum(x[i])`,
#' and `x[i]` is the abundance of species i within a given site.
#' For sites without any species present (Total abundance == 1), Simpson's diversity
#' as calculated here will be undefined (`NaN`).
#' Internal function for use in `hib_calc_div_within()`
#'
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within `hib_calc_div_within()`.
#' @param msg `logical`, indicating if default progress message should be 
#' displayed or not. For use with calculation of Hill-Simpson diversity.
#'
#' @return tibble
#'
#' @references 
#' Simpson, E. H. (1949). Measurement of Diversity. 
#' Nature, 163(4148), 688–688. https://doi.org/10.1038/163688a0
#'
#' @keywords database, processing, diversity
#'
#' @seealso
#' [hib_calc_div_within()], `hib_calc_pie()`
#'
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_simpson(df, "Measurement")
#' } 
#'
#' @noRd
hib_calc_simpson <- function(spp_df, measure_var, msg = TRUE){

  # Assertion checks
  assertthat::assert_that(
    is.data.frame(spp_df),
    is.character(measure_var),
    is.logical(msg)
  )

  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"),
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  # optional message - drop if used for Hill-Simpson
  if (msg){
    message(crayon::white("\tSimpson diversity ...\n"))
  }
  
  # Filter to useable records
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance") 

  simp_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::summarise(Simpson_diversity = 1 - sum( ( !!dplyr::sym(measure_var) / sum(!!dplyr::sym(measure_var)) )^2 ) ) |>
    dplyr::ungroup()
  
  # Ensure all values between 0 and 1!
  simp_df$Simpson_diversity[simp_df$Simpson_diversity < 0] <- 0
  simp_df$Simpson_diversity[simp_df$Simpson_diversity > 1] <- 1
  
  # Above is 'Simpson's Index of Diversity'; 0 = no div, 1 = max div
  # Also known as uncorrected PIE and the Gini-Simpson Index

  # Note: If total abundance == 0, Simpson will be undefined
  # If only 1 species, index value is 0 - i.e no diversity if all one species...
  # NaN for empty community seems acceptable 

  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  out_tib <- dplyr::tibble(Metric = "Simpson_diversity",
                           Data = list(simp_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa) 

  return(out_tib)
}


#' hib_calc_chao1 
#' 
#' @description 
#' Calculate Chao-1 species richness (or diversity) per site. 
#' 
#' @details 
#' Chao-1 is a measure of species richness that also accounts for the prevelance 
#' of rare species (singletons and doubletons).
#' Here, we calculate bias-corrected Chao-1 as: S + (((N_1 * (N_1 - 1)) / (2 * (N_2 + 1)))),
#' where S is species richness, N_1 the number of singleton species, and N_2
#' the number of doubleton species.
#' Internal function for use in `hib_calc_div_within()`.
#' Note, only sites reporting measurements that are suitable for Chao-1 estimation
#' are included. Requirements include integer abundance values.
#' Inspired by code in <https://github.com/timnewbold/predicts-demo/blob/master/predictsFunctions/R/SiteMetrics.R>.
#'
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within `hib_calc_div_within()`.
#' 
#' @return tibble
#' 
#' @references 
#' * Chao, A. (1987). Estimating the Population Size for Capture-Recapture Data 
#' with Unequal Catchability. Biometrics, 43(4), 783–791. https://doi.org/10.2307/2531532
#' * https://github.com/timnewbold/predicts-demo/blob/master/predictsFunctions/R/SiteMetrics.R
#' 
#' @keywords database, processing, diversity
#' 
#' @seealso
#' [hib_calc_div_within()], `hib_calc_sr()`
#' 
#' @examples 
#' \dontrun{
#'   out <- hib_calc_chao1(spp_df, measure_var)
#' } 
#' 
#' @noRd
hib_calc_chao1 <- function(spp_df, measure_var){
  # Assertion checks
  assertthat::assert_that(
    is.data.frame(spp_df),
    is.character(measure_var)
  )
  
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"),
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Diversity_metric_is_suitable_for_Chao"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  message(crayon::white("\tChao1 species richness ...\n"))
  
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance" &
                    Diversity_metric_is_suitable_for_Chao) |>
    # Could this group_by SSBS - if some sites within a study are suitable, why exclude all
    dplyr::group_by(SS) |>
    dplyr::filter(all(floor(!!dplyr::sym(measure_var)) == !!dplyr::sym(measure_var))) |>
    dplyr::ungroup() 

  chao1_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::summarise(Chao1_diversity = sum(!!dplyr::sym(measure_var) > 0) + (((sum(!!dplyr::sym(measure_var) == 1) * (sum(!!dplyr::sym(measure_var) == 1) - 1)) / (2 * (sum(!!dplyr::sym(measure_var) == 2) + 1)))) ) |>
    dplyr::ungroup()

  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  out_tib <- dplyr::tibble(Metric = "Chao1_diversity",
                           Data = list(chao1_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa) 
  
  return(out_tib)
}


#' hib_calc_rsr
#' 
#' @description 
#' Calculate rarefied species richness (RSR) per site.
#' 
#' @details 
#' Rarefied species richness is a measure of species richness that also accounts 
#' for the extent of sampling.
#' Here, we calculate rarefied richness by sampling N individuals from each site,
#' such that N equals the minimum total abundance of a site within a study. All 
#' sites within a study are thus rarefied to equal effort in terms of individuals.
#' Species richness is calculated per site based on the rarefied sample.
#' We repeat the sampling 100 times and take the mean of rarefied species richness 
#' per site.
#' Internal function for use in `hib_calc_div_within()`.
#' Note, only sites reporting measurements that are suitable for Chao-1 estimation
#' are included. Requirements include integer abundance values.
#' Inspired by code in <https://github.com/timnewbold/predicts-demo/blob/master/predictsFunctions/R/SiteMetrics.R>.
#' 
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within `hib_calc_div_within()`.
#' 
#' @return tibble
#'
#' @references 
#' https://github.com/timnewbold/predicts-demo/blob/master/predictsFunctions/R/SiteMetrics.R
#' 
#' @keywords database, processing, diversity
#' 
#' @seealso
#' [hib_calc_div_within()], `hib_calc_sr()`
#' 
#' @examples 
#' \dontrun{
#'   out <- hib_calc_rsr(spp_df, measure_var)
#' } 
#' 
#' @noRd
hib_calc_rsr <- function(spp_df, measure_var){
  # min_N_samp = 10, N_iter = 100

  # TODO: alternative min N to sample (min_N_samp) - 10? or min(c(min_TA, min_N))
  # TODO: alternative (N_iter) - default to 100, but could do 1000? - user specified - store
  
  # Assertion checks
  assertthat::assert_that(
    is.data.frame(spp_df),
    is.character(measure_var)
  )
  
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"),
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Diversity_metric_is_suitable_for_Chao"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  message(crayon::white("\tRarefied species richness ...\n"))

  rarefy_comm <- function(abund_distr, N){
    comm <- rep(as.character(1:length(abund_distr)), abund_distr)
    samp <- sample(comm, N)
    sr <- length(unique(samp))
    return(sr)
  }
  
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance" &
                    Diversity_metric_is_suitable_for_Chao) |>
    dplyr::group_by(SS) |>
    dplyr::filter(all(floor(!!dplyr::sym(measure_var)) == !!dplyr::sym(measure_var))) |>
    dplyr::ungroup()

  rsr_df <- spp_df |>
    dplyr::group_by(SSBS) |>
    dplyr::mutate(TA = sum(!!dplyr::sym(measure_var))) |>
    dplyr::ungroup() |>
    dplyr::group_by(SS) |>
    dplyr::mutate(min_abund = min(TA)) |>
    dplyr::ungroup()
  
  # min_N_samp
  pb <- utils::txtProgressBar(min = 1, max = 100, style = 3)
  
  # N_iter
  rsr_df <- purrr::map_dfr(seq(1,100),
                           function(x, df = rsr_df){
                             tmp <- df |>
                               dplyr::group_by(SSBS, SS) |>
                               dplyr::summarise(sr = rarefy_comm(!!dplyr::sym(measure_var), 
                                                                 min_abund[1]),
                                                min_abund = min_abund[1],
                                                iter = x)
                             
                             utils::setTxtProgressBar(pb, x)
                             tmp
                           }) |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::summarise(Rarefied_species_richness = mean(sr),
                     min_abund = min_abund[1])

  close(pb)
  
  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  out_tib <- dplyr::tibble(Metric = "Rarefied_species_richness",
                           Data = list(rsr_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa)
  
  return(out_tib)

}



#' hib_calc_ma
#' 
#' @description 
#' Calculate arithmetic mean abundance for each site measuring abundance.
#' 
#' @details 
#' Mean abundance = N/S, where N is the total abundance of species within 
#' a given site, and S is the number of species present.
#' Rescaled mean abundance = `Mean abundance[j,k]`/max(`Mean abundance[j,k]`), where
#' the mean abundance in eack site j is divided by the maximum mean abundance 
#' in study k.
#' Internal function for use in `hib_calc_div_within()`.
#' 
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within `hib_calc_div_within()`.
#' 
#' @return tibble
#' 
#' @keywords database, processing, diversity
#' 
#' @seealso
#' [hib_calc_div_within()]
#' 
#' @examples 
#' \dontrun{
#'   out <- hib_calc_ma(spp_df, mv)
#' }
#'  
#' @noRd
hib_calc_ma <- function(spp_df, measure_var){
  # Assertion checks
  assertthat::assert_that(
    is.data.frame(spp_df),
    is.character(measure_var)
  )

  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), 
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  message(crayon::white("\tMean abundance ...\n"))

  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance")

  abund_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::summarise(Mean_abundance = dplyr::case_when(all(!!dplyr::sym(measure_var) == 0) ~ 0,
                                                       T ~ sum(!!dplyr::sym(measure_var)) / sum(!!dplyr::sym(measure_var) > 0) ) ) |>
    dplyr::ungroup() |>
    # Re-scale abundance to 0-1 per study (SS) 
    dplyr::group_by(SS) |>
    # if max Mean_abundance is 0..., leave as is..
    dplyr::mutate(max_ma = max(Mean_abundance),
                  Rescaled_mean_abundance = dplyr::case_when(max_ma == 0 ~ Mean_abundance,
                                                              T ~ Mean_abundance / max_ma)) |>
    dplyr::ungroup() |>
    dplyr::select(-max_ma)

  # Ensure all rescaled values between 0 and 1
  abund_df$Rescaled_mean_abundance[abund_df$Rescaled_mean_abundance < 0] <- 0
  abund_df$Rescaled_mean_abundance[abund_df$Rescaled_mean_abundance > 1] <- 1
  
  # Get N species used
  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  out_tib <- dplyr::tibble(Metric = "Mean_abundance",
                           Data = list(abund_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa)
  
  return(out_tib)
}


#' hib_calc_gma
#' 
#' @description 
#' Calculate geometric mean abundance for each site measuring abundance.
#' 
#' @details 
#' Geometric mean abundance (GMA) is calculated as: exp(mean(log(x_i))), where x_i is 
#' the abundance of species i in a given site.
#' Species with 0 abundance are excluded from calculations.
#' We assign empty sites (0 species richness) a GMA value of 0.
#' Rescaled geometric mean abundance = `Geometric mean abundance[j,k]`/max(`Geometric mean abundance[j,k]`), where
#' the geometric mean abundance in eack site j is divided by the maximum geometric mean abundance 
#' in study k.
#' Internal function for use in `hib_calc_div_within()`.
#' 
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within `hib_calc_div_within()`.
#' 
#' @return tibble
#' 
#' @keywords database, processing
#' 
#' @seealso
#' [hib_calc_div_within()]
#' 
#' @examples 
#' \dontrun{
#'   out <- hib_calc_gma(spp_df, measure_var)
#' }
#'  
#' @noRd
hib_calc_gma <- function(spp_df, measure_var){
  
  # Assertion checks
  assertthat::assert_that(
    is.data.frame(spp_df),
    is.character(measure_var)
  )
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), # obj$species
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  message(crayon::white("\tGeometric mean abundance ...\n"))
  
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance")

  abund_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::filter(!!dplyr::sym(measure_var) > 0) |>
    dplyr::summarise(Geometric_mean_abundance = exp(mean(log(!!dplyr::sym(measure_var)))) ) |>
    dplyr::ungroup()
  
  zero_abund_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::summarise(TA = sum(!!dplyr::sym(measure_var))) |>
    dplyr::ungroup() |>
    dplyr::filter(TA == 0) |>
    dplyr::rename(Geometric_mean_abundance = TA)
  
  abund_df <- abund_df |>
    dplyr::bind_rows(zero_abund_df) |>
    # Re-scale abundance to 0-1 per study (SS) 
    dplyr::group_by(SS) |>
    # if max Geometric_mean_abundance is 0..., leave as is..
    dplyr::mutate(max_gma = max(Geometric_mean_abundance),
                  Rescaled_geometric_mean_abundance = dplyr::case_when(max_gma == 0 ~ Geometric_mean_abundance,
                                                                       T ~ Geometric_mean_abundance / max_gma)) |>
    dplyr::ungroup() |>
    dplyr::select(-max_gma)
  
  # Ensure all rescaled values between 0 and 1
  abund_df$Rescaled_geometric_mean_abundance[abund_df$Rescaled_geometric_mean_abundance < 0] <- 0
  abund_df$Rescaled_geometric_mean_abundance[abund_df$Rescaled_geometric_mean_abundance > 1] <- 1
  
  # Get N species used
  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  out_tib <- dplyr::tibble(Metric = "Geometric_mean_abundance",
                           Data = list(abund_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa)
  
  return(out_tib)
}


#' hib_calc_shannon
#' 
#' @description 
#' Calculate Shannon diversity for each site measuring abundance.
#' 
#' @details 
#' Shannon diversity is calculated as: -1 * sum( p_i * log(p_i)), where p_i is 
#' the relative abundance of species i in the focal site.
#' Empty sites have undefined Shannon diversity (`NA`).
#' Internal function for use in `hib_calc_div_within()`.
#' 
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within `hib_calc_div_within()`.
#' @param msg Logical, indicating if default progress message should be displayed 
#' or not.
#' 
#' @return tibble
#' 
#' @references 
#' Shannon, C. E. (1948). A mathematical theory of communication. The Bell System 
#' Technical Journal, 27(3), 379–423. 
#' 
#' @keywords database, processing, diversity
#' 
#' @examples 
#' \dontrun{
#'   out <- hib_calc_shannon(spp_df, mv)
#' } 
#' 
#' @noRd
hib_calc_shannon <- function(spp_df, measure_var, msg = TRUE){
  # Assertion checks
  assertthat::assert_that(
    is.data.frame(spp_df),
    is.character(measure_var),
    is.logical(msg)
  )
  
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), 
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )

  if (msg){
    message(crayon::white("\tShannon diversity ...\n"))  
  }
  
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance")

  div_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::mutate(p_i = !!dplyr::sym(measure_var) / sum(!!dplyr::sym(measure_var)) ) |>
    # If all p_i == 0, or NA, Shannon_diversity = NA
    dplyr::summarise(Shannon_diversity =  dplyr::case_when(all(p_i == 0) ~ NA,
                                                          all(is.na(p_i)) ~ NA,
                                                           T ~ -1 * sum( p_i * log(p_i), na.rm = TRUE ) ),
                    logN = log(sum(!!dplyr::sym(measure_var) > 0)) ) |>
    dplyr::ungroup()

  # Get N species used
  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  out_tib <- dplyr::tibble(Metric = "Shannon_diversity",
                           Data = list(div_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa)
  
  return(out_tib)
}


#' hib_calc_pie
#' 
#' @description 
#' Calculate PIE (Hulbert's (1971) Probability of Interspecific Encounter), and
#' its species equivalent (sPIE).
#' 
#' @details 
#' PIE is calculated as: (N/(N-1)) * (1 - sum(p_i^2)), where: p_i is the relative
#' abundance of species i in the focal site, and N is the total abundance of organisms
#' in the focal site.
#' This measurement is equivalent to Simpson's Index but includes a correction
#' factor base on the total abundance of organisms.
#' PIE will be undefined (`NA`) if the total abundance of a site is <1 (e.g., non
#' count data), if only a single individual is present, if only singleton species
#' are present.
#' sPIE is calculated as: (1 / (1 - PIE)), equivalent to Hill-Simpson diversity.
#' sPIE is undefined (`NA`) when PIE is undefined, and when PIE == 1.
#' Internal function for use in `hib_calc_div_within()`.
#' 
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within `hib_calc_div_within()`.
#' 
#' @return tibble
#' 
#' @references 
#' Hurlbert, S. H. (1971). The Nonconcept of Species Diversity: A Critique and 
#' Alternative Parameters. Ecology, 52(4), 577–586. https://doi.org/10.2307/1934145
#' 
#' @keywords database, processing, diversity
#' 
#' @seealso
#' [hib_calc_div_within()], `hib_calc_simpson()`
#' 
#' @examples 
#' \dontrun{
#'   out <- hib_calc_pie(df, mv)
#' }
#'  
#' @noRd
hib_calc_pie <- function(spp_df, measure_var){
  # Assertion checks
  assertthat::assert_that(
    is.data.frame(spp_df),
    is.character(measure_var)
  )
  
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), 
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  message(crayon::white("\tProbability of Interspecific Encounter (PIE) ...\n"))
  
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance") 
  
  # If non-counts in abundance - warning / flag in output
  if (any(spp_df[,measure_var] %% 1 != 0, na.rm = TRUE)){
    message(crayon::yellow("Warning: use of non-count abudance data could lead to issues calculating PIE..."))
  }
  # if measure_var != "Measurement" - print a warning and change measure_var

  # if (any(grepl("percent cover"), spp_df$Diversity_metric))
  
  
  # TO CHECK - can PIE be calculated on effort corrected measures? -
  # Mobr code suggests not...
  # Should really be counts 
  # Also only really works for counts...
  # TODO: fixes for empty comms - currently, if all 0, PIE == NAN, should it be 0?
  # TODO: fixes for comm of singletons, currently PIE == 1, sPIE == Inf, zero?, 1 - undefined
  
  # NOTES: in other functions, if N (total abund) < 2, PIE == NA
  # If PIE == 1, sPIE == NA
  
  # In MOBR, if 0/all zero, PIE == 0, sPIE == 0
  # In MOBR, if only 1 indiv, PIE == NA, sPIE == NA
  # In MOBR, if multiple singletons, PIE == 1, sPIE == NA
  
  # simpson ~ PIE
  # In vegan, if simpson, and all 0, output == 1, thus inv_simpson == Inf!!
  # In vegan, if simpson, and single 1, output == 0, thus inv_simpson == 1!!
  # Vegan treats NAs as 0s, ish, drops NAs from summation...
  # How to address 
  
  pie_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::summarise(adj_val = sum( !!dplyr::sym(measure_var) ) / ( sum( !!dplyr::sym(measure_var) ) - 1 ),
                     PIE_ = 1 - sum( ( !!dplyr::sym(measure_var) / sum(!!dplyr::sym(measure_var)) )^2 ),
                     PIE = adj_val * PIE_,
                     sPIE = (1 / (1 - PIE)),
                     
                     N_spp = sum(!!dplyr::sym(measure_var) > 0 ),
                     Tot_abund =  sum(!!dplyr::sym(measure_var)) ) |>
    dplyr::ungroup() 
  
  # Ensure all PIE values between 0 and 1
  pie_df$PIE[pie_df$PIE < 0] <- 0
  pie_df$PIE[pie_df$PIE > 1] <- 1
  
  # Based on MOBR, do the following
  # If N (total abund) < 2, PIE should be NA ?, or 0
  # Tot_abund, N_spp
  pie_df$PIE[pie_df$Tot_abund <= 1] <- NA
  pie_df$PIE[pie_df$Tot_abund == 0] <- 0
  
  if (any(pie_df$Tot_abund == 1)){
    message(crayon::yellow("Warning: some sites contain only one individual, PIE will be undefined (NA) ..."))
  }
  if ( any( ((pie_df$Tot_abund <= 1)  &  (pie_df$Tot_abund > 0)) ) ){
    message(crayon::yellow("Warning: some sites contain fewer than one individual, PIE will be undefined (NA) ..."))
  }
  
  if ( any( sapply(pie_df$PIE, function(x){isTRUE(all.equal(x, 1))}) ) ){
    message(crayon::yellow("Warning: some sites contain only singletons, sPIE will be undefined (NA) ..."))
  }
  
  # If PIE == 1, sPIE should be NA  - 
  # If PIE == 0, sPIE should be 0
  pie_df$sPIE[is.na(pie_df$PIE)] <- NA
  pie_df$sPIE[pie_df$PIE == 1] <- NA
  pie_df$sPIE[pie_df$PIE == 0] <- 0
  
  pie_df$sPIE[sapply(pie_df$PIE, function(x){isTRUE(all.equal(x, 1))})] <- NA
  pie_df$sPIE[sapply(pie_df$PIE, function(x){isTRUE(all.equal(x, 0))})] <- 0
  
  # pie_df$sPIE[is.infinite(pie_df$sPIE)] <- NA
  
  pie_df <- pie_df |> 
    dplyr::select(SSBS, SS, PIE, sPIE)
  
  # Get N species used
  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  out_tib <- dplyr::tibble(Metric = "PIE",
                           Data = list(pie_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa)
  
  return(out_tib)
}


#' get_eq_ev
#'
#' @description 
#' Calculcate equalness (and evenness) of a site.
#' 
#' @details 
#' Used internally to calculate equalness and/or evenness as in Gregorius & 
#' Gillet, 2024.
#' Here, restricted to returning equalness. Called by `hib_calc_equalness()`.
#' 
#' @param f Numeric vector with absolute frequencies
#' @param q Order of Hill numbers
#'
#' @return numeric
#'
#' @references H-R Gregorius, EM Gillet. 2024.
#' On the diversity-based measures of equalness and evenness.
#' Methods in Ecology and Evolution,
#' https://doi.org/10.1111/2041-210X.14265
#' 
#' @seealso
#' `hib_calc_equalness()`, [hib_calc_div_within()]
#'
#' @examples
#' \dontrun{
#'   f = c( 10, 20, 30, 40, 50, 60 )
#'   q = 1
#'   get_eq_ev( f, q )
#' }
#' @noRd
get_eq_ev <- function(f, q = 1){
  assertthat::assert_that(
    is.numeric(f),
    is.numeric(q)
  )

  # number of types with a positive frequency
  n = length( f[f > 0] )
  
  # decreasing relative frequencies
  f = f / sum( f )
  f = sort(f, decreasing = TRUE)
  
  # relative step heights
  h = ( f - c( f[-1], 0 ) ) / f[1]
  h = h / sum( h )
  
  # diversity for frequencies (vf) and for step heights (vh)
  # one single type
  vf = 1.0 # diversity for frequency distribution
  vh = 1.0 # diversity for step-height distribution
  if ( n > 1 ){
    if ( q == 1 ){
      # effective number of Shannon entropy (Hill number of order 1)
      vf = exp( -sum( f * log( f ) , na.rm = TRUE ) ) # Need to ignore 0 in f...
      hp = h[h > 0] # only different heights contribute diversity
      vh = exp( -sum( hp * log( hp ) ) )
    } else{
      # other Hill numbers
      vf = (sum(f^q))^(1./(1.-q))
      vh = (sum(h^q))^(1./(1.-q))
    }
  }
  
  # evenness
  ef = 0.0 # measure of evenness
  eh = 1.0 # measure of equalness
  
  if ( n > 1 ){
    ef = (vf - 1.0) / (n - 1.0)
    eh = (n - vh) / (n - 1.0)

    if (ef < 0.0001){
      ef = 0
    }
    if (eh < 0.0001){
      eh = 0
    }
  }

  out <- eh
  return(out)
}


#' hib_calc_equalness 
#' 
#' @description 
#' Calculate equalness per site, following Gregorius and Gillet, 2024.
#' 
#' @details 
#' Equalness is defined as: (S - D) / (S - 1), where S is species richness, and
#' D is Hill-Shannon (exp(Shannon)).
#' Empty communities, single species communities, and even communities (all species
#' have the same abundance) have an equalness of 1.
#' 
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within `hib_calc_div_within()`.
#' 
#' @return tibble
#' 
#' @references 
#' H-R Gregorius, EM Gillet. 2024. On the diversity-based measures of equalness 
#' and evenness. Methods in Ecology and Evolution, https://doi.org/10.1111/2041-210X.14265
#' 
#' @keywords database, processing
#' 
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_equalness(df, "Measurement")
#' }
#'  
#' @noRd
hib_calc_equalness <- function(spp_df, measure_var){
  
  # NOTE: could use evenness from Gregorious? - similar to default here, but - 1 in 
  # numerator and denominator
  
  # Assertion checks
  assertthat::assert_that(
    is.data.frame(spp_df),
    is.character(measure_var)
  )
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), # obj$species
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  message(crayon::white("\tEqualness ...\n"))

  # Filter to useable records
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance") 
  
  eq_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::summarise(Equalness = get_eq_ev( !!dplyr::sym(measure_var)) ) |>
    dplyr::ungroup()
  
  spp_info <- get_n_spp_taxa(spp_df, measure_var)
  
  out_tib <- dplyr::tibble(Metric = "Equalness",
                           Data = list(eq_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa) 
  
  return(out_tib)
}


#' Calculate within-site diversity metrics
#' 
#' @description 
#' This function calculates multiple, user specified measures of site-level 
#' (taxonomic) diversity for sites in the associated HibDatabase.
#' 
#' @details 
#' As is standard for R6 Classes, the HibDatabase is modified in-place.
#' If effort correction (`hib_effort_correction()`) has previously been conducted, 
#' diversity metrics will be estimated using effort corrected measurements.
#' 
#' @param obj A [`HibDatabase`] object.
#' @param metrics Character vector, any/all of; "Total_abundance", "Species_richness", 
#' "Chao1_diversity",  "Simpson_diversity", "Rarefied_species_richness",
#' "Mean_abundance", "Geometric_mean_abundance", "Shannon_diversity",
#' "HillShannon_diversity", "HillSimpson_diversity", "PIE", "Equalness". 
#' Metrics can be provided in upper or lower case and as unambiguous abbreviations 
#' (internal processing is done using lower-case).
#' 
#' @return HibDatabase with updated '$diversity_within' tibble
#' 
#' @keywords database, processing
#' 
#' @examples 
#' \dontrun{
#'   hib_db <- hib_calc_div_within(hib_db, c("Total_abundance", "Species_richness"))
#' }
#'  
#' @export
hib_calc_div_within <- function(obj, metrics){
  
  # Assertions
  assertthat::assert_that(
    inherits(obj, "HibDatabase"),
    assertthat::has_name(obj$species, "SSBS"),
    assertthat::has_name(obj$species, "SS"),
    assertthat::has_name(obj$species, "Measurement") || assertthat::has_name(obj$species, "Effort_corrected_measurement"),
    assertthat::has_name(obj$species, "Diversity_metric_type"),
    assertthat::has_name(obj$species, "Best_guess_binomial"),
    is.character(metrics)
  )

  # TODO:
  # Custom metric - 
  # use "custom_metric" in metrics vec
  # add a custom_metric_spec arg - list/named vector...
  # have custom_metric_func, custom_metric_name 
  
  
  metrics <- tolower(metrics)
  # Match supplied metrics to options.  
  metrics <- match.arg(metrics,
                       c("total_abundance", "species_richness", "chao1_diversity",  
                         "simpson_diversity", "rarefied_species_richness", #),
                         "mean_abundance", "geometric_mean_abundance",
                         "shannon_diversity", "evenness",
                         "hillshannon_diversity", "hillsimpson_diversity",
                         "pie", "equalness"),
                       several.ok = T)
  metrics <- unique(metrics)

  # %%%
  # if metric is already in hib$diversity_within,
  if (!is.null(obj$diversity_within)){
    metrics_tc <- tools::toTitleCase(metrics)
    
    # %%%
    # "HillShannon_diversity", "HillSimpson_diversity", "PIE" need to be
    # separately edited from lowercase to 'tc' equivalents
    metrics_tc <- metrics
    metrics_tc[metrics_tc == "hillshannon_diversity"] <- "HillShannon_diversity"
    metrics_tc[metrics_tc == "hillsimpson_diversity"] <- "HillSimpson_diversity"
    metrics_tc[metrics_tc == "pie"] <- "PIE"
    metrics_tc <- tools::toTitleCase(metrics_tc)
    # %%%
    new_metrics_idx <- which(!(metrics_tc %in% obj$diversity_within$Metric))
    # drop from metrics and print the dropped metrics...
    metrics <- metrics[new_metrics_idx]
    dropped_metrics <- metrics_tc[metrics_tc %in% obj$diversity_within$Metric]
    if (length(dropped_metrics) > 0){
      message(crayon::yellow(paste0("\nWarning: the following metrics have already been ",
                                  "added to the HibDatabase and so will not be ",
                                  "recalculated;\n\t",
                                  paste(dropped_metrics, collapse = ", ")
                                  # "\n"
                                  )))
    }
    if (length(metrics) == 0){
        stop("No new biodiversity metrics to calculate.", call. = F)
    }
  }
  # %%%

  message(crayon::green("\nCalculating within-site biodiversity metrics ...\n"))
  
  if (obj$processing_flags$measures_effort_corrected == T){
    measure_var <- "Effort_corrected_measurement"
  } 
  if (obj$processing_flags$measures_effort_corrected == F){
    measure_var <- "Measurement"
  }
  # Check for NAs in Measurement/Effor_corrected_measurement...
  # Set to 0 if any have crept in...
  if (sum(is.na(obj$species[,measure_var])) > 0){
    message(crayon::yellow(paste0("\nWarning: NAs detected in '", measure_var, "'.  Setting NAs to zeros.")))
    obj$species[is.na(obj$species[,measure_var]), measure_var] <- 0
  }

  assertthat::assert_that(
    sum(is.na(obj$species[,measure_var])) == 0,
    is.numeric(obj$species[[measure_var]])
  )

  out_tib <- purrr::map(.x = metrics, 
                       .f = function(metric = .x, db_obj = obj,
                                    spp_df = obj$species, mv = measure_var){
                         if (metric == "total_abundance"){
                           return(hib_calc_ta(spp_df, mv))
                         }
                         if (metric == "species_richness"){
                           return(hib_calc_sr(spp_df, mv))
                         }
                         if (metric == "chao1_diversity"){
                           return(hib_calc_chao1(spp_df, mv))
                         }
                         if (metric == "simpson_diversity"){
                           return(hib_calc_simpson(spp_df, mv))
                         }
                         if (metric == "rarefied_species_richness"){
                           return(hib_calc_rsr(spp_df, mv))
                         }
                         if (metric == "mean_abundance"){
                           return(hib_calc_ma(spp_df, mv))
                         }
                         if (metric == "geometric_mean_abundance"){
                           return(hib_calc_gma(spp_df, mv))
                         }
                         
                         if (metric == "shannon_diversity"){
                           tmp <- hib_calc_shannon(spp_df, mv)
                           tmp$Data[[1]] <- tmp$Data[[1]] |>
                             dplyr::select(SSBS, SS, Shannon_diversity)
                           return(tmp)
                         }
                         # %%%
                         if (metric == "pie"){
                           return(hib_calc_pie(spp_df, mv))
                         }
                         # %%%
                         # Evenness (Shannon/log(N_spp))
                         if (metric == "evenness"){
                           message(crayon::white("\tEvenness ...\n"))

                           tmp <- hib_calc_shannon(spp_df, mv, msg = FALSE)
                           tmp$Metric <- "Evenness"
                           tmp$Data[[1]]$Evenness <- tmp$Data[[1]]$Shannon_diversity / tmp$Data[[1]]$logN
                           tmp$Data[[1]] <- tmp$Data[[1]] |>
                             dplyr::select(SSBS, SS, Evenness, logN) |>
                             # If only 1 species, set evenness to 0 to avoid division by 0 error 
                             dplyr::mutate(Evenness = dplyr::case_when(logN == 0 ~ 0, 
                                                                       T ~ Evenness)) |> 
                             # Ensure Evenness is 0-1
                             dplyr::mutate(Evenness = dplyr::case_when(Evenness < 0 ~ 0,
                                                                       Evenness > 1 ~ 1,
                                                                       T ~ Evenness))
                           
                           return(tmp)
                         }
                         # TODO:
                         # Alternative evenness? - (1/Simpson)/SR
                         # Coverage-based richness - N spp if 80% of comm sampled
                         
                         # %%%
                         # Hill-Shannon: exp(shannon)
                         if (metric == "hillshannon_diversity"){
                           message(crayon::white("\tHill-Shannon diversity ...\n"))

                           tmp <- hib_calc_shannon(spp_df, mv, msg = FALSE)
                           tmp$Metric <- "HillShannon_diversity"
                           tmp$Data[[1]]$HillShannon_diversity <- exp(tmp$Data[[1]]$Shannon_diversity)
                           
                           tmp$Data[[1]] <- tmp$Data[[1]] |>
                             dplyr::select(SSBS, SS, HillShannon_diversity)
                           return(tmp)
                         }
                         
                         # Hill-Simpson: 1/simpson 
                         if (metric == "hillsimpson_diversity"){
                           message(crayon::white("\tHill-Simpson diversity ...\n"))

                           tmp <- hib_calc_simpson(spp_df, mv, msg = FALSE)
                           tmp$Metric <- "HillSimpson_diversity"
                           tmp$Data[[1]]$HillSimpson_diversity <- 1 / ( 1 - tmp$Data[[1]]$Simpson_diversity )
                           
                           tmp$Data[[1]] <- tmp$Data[[1]] |>
                             dplyr::select(SSBS, SS, HillSimpson_diversity)
                           return(tmp)
                         }
                         if (metric == "equalness"){
                           return(hib_calc_equalness(spp_df, mv))
                         }
                         # %%%
                         }) |> dplyr::bind_rows()
  
  # Effort correction flag
  out_tib$Flag <- paste0("Effort corrected: ", obj$processing_flags$measures_effort_corrected)

  div_empty <- is.null(obj$diversity_within)
  
  if (!div_empty){
    obj$diversity_within <- dplyr::bind_rows(obj$diversity_within, out_tib)
  }  
  if (div_empty){
    obj$diversity_within <- out_tib
  }
  
  obj$last_update <- lubridate::now()

  # checkself?

  return(obj)
}




#' hib_calc_veg
#' 
#' @description 
#' This function calculates vegan-based dissimilarity metrics between sites.
#' 
#' @details 
#' Function is largley a wrapper for [`vegan::vegdist()`], with some additional
#' pre-processing of inputs and post-processing of outputs.
#' 
#' @param spp_site_data A data.frame/tibble, with species as columns and sites as
#' rows.
#' @param metric Character specifying metric to calculate - see [`vegan::vegdist()`].
#' @param binary Logical - should distances be based on presence/absence.
#' @param rel_abund Logical - should relative abundance values be used in distance
#' calculations.
#' @param rescale Logical - should distance be rescaled relative to that metric's
#' theoretical maximum (if applicable and not 1).
#' @param div_metric_type A character describing the type of biodiversity records
#' at the site of interest - one of "Occurrence" or "Abundance". Set by `hib_cacl_div_between()`.
#' 
#' @return data.frame
#' 
#' @keywords database, processing, diversity
#' 
#' @seealso
#' [hib_calc_div_between()], [vegan::vegdist()]
#' 
#' @examples 
#' \dontrun{
#'   metric_val <- hib_calc_veg(spp_site_data, "euclidean")
#' }
#'  
#' @noRd
hib_calc_veg <- function(spp_site_data, 
                         metric,
                         binary,
                         rel_abund,
                         rescale,
                         div_metric_type
                         ){
  
  assertthat::assert_that(
    is.data.frame(spp_site_data),
    is.character(metric),
    is.logical(binary) || is.null(binary),
    is.logical(rel_abund) || is.null(rel_abund),
    is.logical(rescale) || is.null(rescale),
    is.character(div_metric_type)
  )
  
  metric <- match.arg(metric, 
                      c("manhattan", "euclidean", "canberra", "clark", "bray", 
                        "kulczynski", "jaccard", "gower", "altGower", 
                        "morisita", "horn", "mountford", "raup", "binomial", 
                        "chao", "cao", "mahalanobis", "chisq", "chord", 
                        "hellinger", "aitchison", "robust.aitchison"), 
                      several.ok = FALSE)
  
  # TODO: NA checks...

  if (rel_abund){
    row_abund <- rowSums(spp_site_data, na.rm = TRUE)
    row_abund[row_abund==0] <- 1
    spp_site_data <- spp_site_data/row_abund
  }
  
  # TODO: which metrics do/don't work with zeros in one/both sites...
  # aitchison

  # if binary is false, and div_metr_type is Occurrence
  # Skip /make an 'empty' dist_df (SSBS1, SSBS2, Distance, metric_nm, Contrast)
  if ((binary == F) & (div_metric_type == "Occurrence")){
    return(NULL)
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
  
  dist_df <- dist_df |>
    dplyr::filter(SSBS1 != SSBS2) |> 
    dplyr::mutate(metric_nm = tools::toTitleCase(metric),
                  contrast = paste(SSBS1, SSBS2, sep = "_"))
  
  if (rescale & !is.na(max_dist)){
    dist_df$Distance <- dist_df$Distance / max_dist
    # Ensure no rounding errors- 0-1
    dist_df$Distance[dist_df$Distance > 1] <- 1
  }

  # Adding rescaling of euclidian if using relative abundance
  if (rescale & rel_abund & metric == "euclidean"){
    dist_df$Distance <- dist_df$Distance / sqrt(2)
    # Ensure no rounding errors- 0-1
    dist_df$Distance[dist_df$Distance > 1] <- 1
  }
  
  # Ensure no rounding errors- > 0 
  dist_df$Distance[dist_df$Distance < 0] <- 0
  
  # Ideally, if rescale is T and max_dist is NA, edit rescale to be F in output
  
  # TODO: if as_sim, convert to similarity score (not possible for unscaled distances...)

  return(dist_df)
}


#' calc_hypervolume
#' 
#' @description 
#' Calculate hypervolume-based dissimilarity/beta-diversity given two communities.
#' 
#' @details 
#' Expects two vectors, one each for the abundance/occurence of species at two sites
#' assuming that values at the same index (e.g., `x[i]`, `y[i]`) correspond to the
#' same species/taxonomic entity.
#' Output is the normalised hypervolume.
#' 
#' Internal function for use in `hib_calc_nonveg()`, called in turn by `hib_calc_div_between()`.
#' 
#' @param spp_site_matrix Numeric matrix containing abundance/occurence of species 
#' at two sites. Sites are rows, species are columns.
#' 
#' @return numeric
#' 
#' @references 
#' Song, Chuliang, Muyang Lu, Joseph R. Bennett, Benjamin Gilbert, Marie-Josée 
#' Fortin, and Andrew Gonzalez (2025). A Geometric Approach to Beta Diversity. 
#' Ecological Monographs 95(1): e70008. https://doi.org/10.1002/ecm.70008
#'
#' @keywords database, processing, diversity
#' 
#' @seealso
#' [hib_calc_div_between()], `hib_calc_nonveg()`
#' 
#' @examples 
#' \dontrun{
#'   hypervolume <- calc_hypervolume(spp_site_matrix)
#' }
#'  
#' @noRd
calc_hypervolume <- function(spp_site_matrix){
  assertthat::assert_that(
    is.data.frame(spp_site_matrix)
  )

  spp_site_matrix <- t(spp_site_matrix)
    
  # d <- min(dim(spp_site_matrix))

  # add origin to data
  spp_site_matrix <- rbind(spp_site_matrix,
                         rep(0, ncol(spp_site_matrix)))

  vol <- tryCatch(
    {
      geometry::convhulln(spp_site_matrix, 
                          output.options = TRUE)$vol
    },
    error = function(e) {
      # TODO: Would NA be a better error value??
      0
    }
  )
  
  # Scale based on min(Nspp, Ncomm) - assume to be 2. 
  # If only 1 species, will return 0 anyway
  vol <- sqrt(vol)

  return(vol)
}


#' calc_cosine_sim 
#' 
#' @description 
#' Calculates Cosine similarity between the (relative) abundances/occurences of 
#' species at two sites.
#' 
#' @details 
#' Expects two vectors, one each for the abundance/occurence of species at two sites
#' assuming that values at the same index (e.g., `x1[i]`, `x2[i]`) correspond to the
#' same species/taxonomic entity.
#' Cosine similarity is calculated as: sum( x1 * x2 ) / ( sqrt( sum( x1^2 )) * sqrt( sum( x2^2 ) ) )
#' Note, that by default, this measure is a similarity, contrasting with most 
#' other between-site measures, which are dissimilarity/distances. Cosine similarity is
#' therefore post-processed to it's inverse ( 1 - Cosine ) in `hib_calc_nonveg()`.
#' Internal function for use in `hib_calc_nonveg()`, called in turn by `hib_calc_div_between()`.
#' 
#' @param x1 Numeric vector for abundance/occurence of species at site 1.
#' @param x2 Numeric vector for abundance/occurence of species at site 2.
#' 
#' @return numeric
#' 
#' @keywords database, processing, diversity
#' 
#' @seealso
#' [hib_calc_div_between()], `hib_calc_nonveg()`
#' 
#' @examples 
#' \dontrun{
#'   cosine_sim <- calc_cosine_sim(c(0,1,1), c(1,0,1))
#' }
#'  
#' @noRd
calc_cosine_sim <- function(x1, x2){
  assertthat::assert_that(
    length(x1) == length(x2),
    is.numeric(x1),
    is.numeric(x2)
  )

  cosine <- sum( x1 * x2 , na.rm = TRUE) / ( sqrt( sum( x1^2 , na.rm = TRUE)) * sqrt( sum( x2^2 , na.rm = TRUE) ) )
  # TODO: to address div by zero if one vector is all 0, could use a small value 
  # "epsilon" (1e-8) to denomintor if zero
  # if one of the sites doesn't have any individuals/species in it
  # i.e. the row sum is 0
  
  # Could set cosine to 0 if one site is empty.
  # And to NA if both sites are?


  # ensure all values between 0 and 1
  if (!is.na(cosine)){
    if (cosine < 0){ cosine <- 0}
    if (cosine > 1){ cosine <- 1}  
  }
  
  return(cosine)
}


#' calc_horn_overlap
#' 
#' @description 
#' Calculate Horn's overlap (1966) metric given two communities. 
#' 
#' @details 
#' Expects two vectors, one each for the abundance/occurence of species at two sites
#' assuming that values at the same index (e.g., `x[i]`, `y[i]`) correspond to the
#' same species/taxonomic entity.
#' Note, that by default, Horn's overlap is a similarity, contrasting with most 
#' other between-site measures, which are dissimilarity/distances. Overlap is 
#' therefore post-processed to it's inverse ( 1 - overlap ) in `hib_calc_nonveg()`.
#' Internal function for use in `hib_calc_nonveg()`, called in turn by `hib_calc_div_between()`.
#' 
#' @param x Numeric vector for abundance/occurence of species at site 1.
#' @param y Numeric vector for abundance/occurence of species at site 2.
#' 
#' @return numeric
#' 
#' @keywords database, processing, diversity
#' 
#' @seealso
#' [hib_calc_div_between()], `hib_calc_nonveg()`
#' 
#' @examples 
#' \dontrun{
#'   horn_overlap <- calc_horn_overlap(c(0,1,1), c(1,0,1))
#' }
#'  
#' @noRd
calc_horn_overlap <- function(x, y){
  assertthat::assert_that(
    length(x) == length(y),
    is.numeric(x),
    is.numeric(y)
  )
  
  num1 <- sum( (x + y) * log(x + y), na.rm = T )
  num2 <- sum(x * log(x), na.rm = T)
  num3 <- sum(y * log(y), na.rm = T)
  
  num <- num1 - num2 - num3
  
  x_sum <- sum(x, na.rm = T)
  y_sum <- sum(y, na.rm = T)
  xy_sum <- x_sum + y_sum
  
  denom1 <- xy_sum * log(xy_sum)
  denom2 <- x_sum * log(x_sum)
  denom3 <- y_sum * log(y_sum)
  
  denom <- denom1 - denom2 - denom3
  
  horn_overlap <- round(num/denom, 5)
  
  if (out < 0){
    out <- 0
  }
  if (out > 1){
    out <- 1
  }
  return(out)
}


#' calc_minkowski
#' 
#' @description 
#' Calculate normalised form sof Minkowski's distance given two communities. In 
#' the default setting this corresponds to a normalised Euclidean distance (alpha = 2).
#' 
#' @details 
#' Expects two vectors, one each for the abundance/occurence of species at two sites
#' assuming that values at the same index (e.g., `x[i]`, `y[i]`) correspond to the
#' same species/taxonomic entity.
#' Output is NUM/DEN, where NUM is calculated as:
#' sum(|x - y|^alpha)^(1/alpha), 
#' and DEN is:
#' sum((x + y)^alpha)^(1/alpha), if method == "D"
#' (sum(x^alpha) + sum(y^alpha))^(1/alpha), if method == "delta"
#' sum(x^alpha)^(1/alpha) + sum(y^alpha)^(1/alpha), if method == "L".
#' 
#' Internal function for use in `hib_calc_nonveg()`, called in turn by `hib_calc_div_between()`.
#' 
#' @param x Numeric vector for abundance/occurence of species at site 1.
#' @param y Numeric vector for abundance/occurence of species at site 2.
#' @param method Character - indicates the denominator for normalisation. One of 
#' "D", "delta", or "L".
#' @param alpha Integer - exponent used in distance calculation.
#' 
#' @return numeric
#' 
#' @references 
#' Ricotta, C., and Pavoine, S. (2022) A new parametric measure of functional 
#' dissimilarity: Bridging the gap between the Bray-Curtis dissimilarity and the 
#' Euclidean distance. Ecological Modelling, Volume 466, April 2022, 109880
#' https://doi.org/10.1016/j.ecolmodel.2022.109880
#'
#' @keywords database, processing, diversity
#' 
#' @seealso
#' [hib_calc_div_between()], `hib_calc_nonveg()`
#' 
#' @examples 
#' \dontrun{
#'   normalised_euclidean <- calc_minkowski(c(0,1,1), c(1,0,1))
#' }
#'  
#' @noRd
calc_minkowski <- function(x, y, method = "D", alpha = 2){
  assertthat::assert_that(
    length(x) == length(y),
    is.numeric(x),
    is.numeric(y),
    is.character(method),
    is.numeric(alpha)
  )

  assertthat::assert_that(
    alpha > 0
  )

  p1 <- sum(abs(x - y)^alpha)^(1/alpha)
  
  if (method == "D"){
    p2 <- sum((x + y)^alpha)^(1/alpha)
  }
  if (method == "delta"){
    p2 <- (sum(x^alpha) + sum(y^alpha))^(1/alpha)
  }
  if (method == "L"){
    p2 <- sum(x^alpha)^(1/alpha) + sum(y^alpha)^(1/alpha)
  }
  
  out <- (p1 / p2)
  
  return(out)
}


#' hib_calc_nonveg
#' 
#' @description 
#' Calculate non-vegan dissimilarity metrics between sites.
#' 
#' @details 
#' Expects a dataframe/tibble with species abundance records at two sites 
#' - sites as rows, species as columns.
#' Largely a wrapper for other functions which calculate the metrics.
#' 
#' @param spp_site_data A data.frame/tibble, with species as columns and sites as
#' rows.
#' @param metric Character specifying metric to calculate - one of "hypervolume", 
#' "cosine", "theta", "baselga_bray_balanced", "baselga_bray_gradient", 
#' "baselga_ruzicka_balanced", "baselga_ruzicka_gradient", "horn_dissim", 
#' "normalised_euclidean".
#' @param binary Logical - should distance be based on presence/absence
#' @param rel_abund Logical - should relative abundance values be used;
#' applicable for hypervolume and baselga_<X> metrics. Cosine and theta are 
#' unaffected by relative abundance v raw abundance.
#' @param rescale Logical - should distance be rescaled relative to that metric's
#' theoretical maximum (if applicable and not 1).
#' @param div_metric_type Character describing the type of biodiversity records
#' at the site of interest - one of "Occurrence" or "Abundance". Set by `hib_cacl_div_between()`.
#' 
#' @keywords database, processing, diversity
#' 
#' @return data.frame
#' 
#' @seealso
#' [hib_calc_div_between()], `hib_calc_veg()`
#' 
#' @examples 
#' \dontrun{
#'   metric_val <- hib_calc_nonveg(spp_site_data, "cosine")
#' }
#'  
#' @noRd
hib_calc_nonveg <- function(spp_site_data, 
                            metric,
                            binary,
                            rel_abund,
                            rescale,
                            div_metric_type
                            ){
  
  assertthat::assert_that(
    is.character(metric),
    is.logical(binary) || is.null(binary),
    is.logical(rel_abund) || is.null(rel_abund),
    is.logical(rescale) || is.null(rescale),
    is.character(div_metric_type)
  )
  
  metric <- match.arg(metric, 
                      c("hypervolume", "cosine", "theta",
                        "baselga_bray_balanced", "baselga_bray_gradient", 
                        "baselga_ruzicka_balanced", "baselga_ruzicka_gradient",
                        "horn_dissim", "normalised_euclidean"),
                      # "baselga_bray","baselga_ruzicka",
                      # "baselga_jaccard", "baselga_sorensen") 
                      several.ok = FALSE)
  # TODO: challenge of single spp communities?? - this would affect metric calcs..??

  # Check for binary
  if ((binary == F) & (div_metric_type == "Occurrence")){
    return(NULL)
  }
  # Convert to relative abundance if needed
  if (rel_abund){
    if (metric %in% c("hypervolume", 
                        "baselga_bray_balanced", "baselga_bray_gradient", 
                        "baselga_ruzicka_balanced", "baselga_ruzicka_gradient",
                        "horn_dissim", "normalised_euclidean")){
      row_abund <- rowSums(spp_site_data, na.rm = TRUE)
      row_abund[row_abund == 0] <- 1
      spp_site_data <- spp_site_data/row_abund
    }
  }
  # Use species occurence if binary specified
  if (binary){
    spp_site_data[spp_site_data > 0] <- 1
  }
  # if!binary, abundance data only is filtered later

  # Calculate hypervolume
  if (metric == "hypervolume"){
    # Should really always use relative abundance! (or binary) - raw abundance doesn't make sense
    # TODO: add a warning if binary == F and rel_abund == F

    # # In Song et al, vol is further rescaled by min(Nspp, Ncomm) - not done here
    # # We simply return normalised hypervolume - implicitly has max of 1. 
    # What about max for binary inputs?? 

    vol <- calc_hypervolume(spp_site_data)

    out <- data.frame("Hypervolume" = vol)
  }
  
  # Calculate cosine/theta
  else if (metric %in% c("cosine", "theta")){
    
    x <- as.matrix(spp_site_data)
    # N rows == 2
    x1 <- x[1,]
    x2 <- x[2,]

    cosine_sim <- calc_cosine_sim(x1, x2)

    if (metric == "cosine"){
      out <- data.frame("Cosine" = 1 - cosine_sim)
    }
    if (metric == "theta"){
      theta <- acos(cosine_sim)
      if (!is.na(theta)){
        if (theta < 0){ 
          theta <- 0
        }  
      }
      
      if (rescale){
        theta <- theta / (pi / 2)
        # Ensure all rescaled values between 0 and 1!
        if (!is.na(theta)){
          if (theta > 1){ 
            theta <- 1
          }  
        }
      }
      # Else ensure 0 - pi/2
      if (!is.na(theta)){
        if (theta > (pi/2)){ 
          theta <- (pi/2)
        }
      }
      out <- data.frame("Theta" = theta)
    }
  }

  # Calculate horn
  else if (metric == "horn_dissim"){
    x <- as.matrix(spp_site_data)
    # N rows == 2
    x1 <- x[1,]
    x2 <- x[2,]
    horn_overlap <- calc_horn_overlap(x1, x2)
    out <- data.frame("Horn_dissim" = 1 - horn_overlap)
  }
  
  # Calculate baselga_X metrics
  else if (metric %in% c("baselga_bray_balanced", "baselga_bray_gradient",
                         "baselga_ruzicka_balanced", "baselga_ruzicka_gradient")){
    
    if (metric %in% c("baselga_bray_balanced", "baselga_bray_gradient")){
      out <- unlist(betapart::beta.pair.abund(spp_site_data, 
                                              index.family = "bray"))
    }
    if (metric %in% c("baselga_ruzicka_balanced", "baselga_ruzicka_gradient")){
      out <- unlist(betapart::beta.pair.abund(spp_site_data, 
                                              index.family = "ruzicka"))
    }
    
    out <- as.data.frame.list(out)

   
    # Rename columns
    if (metric %in% c("baselga_bray_balanced", "baselga_bray_gradient")){
      out <- out |>
        dplyr::rename(baselga_bray_balanced = beta.bray.bal,
                      baselga_bray_gradient = beta.bray.gra)
    }
    if (metric %in% c("baselga_ruzicka_balanced", "baselga_ruzicka_gradient")){
      out <- out |>
        dplyr::rename(baselga_ruzicka_balanced = beta.ruz.bal,
                      baselga_ruzicka_gradient = beta.ruz.gra)
    }
    
    out <- out |>
      dplyr::select(!!metric)
    
    colnames(out) <- tools::toTitleCase(colnames(out))
  }
  
  # normalised_euclidean
  else if (metric == "normalised_euclidean"){
    x <- as.matrix(spp_site_data)
    # N rows == 2
    x1 <- x[1,]
    x2 <- x[2,]
    
    norm_euc <- calc_minkowski(x1, x2)
    out <- data.frame("Normalised_euclidean" = norm_euc)
  }
  # TODO: if as_sim, convert to similarity score (not possible for unscaled distances...)
  
  return(out)
}


#' Calculate between-site (within study) diversity metrics
#' 
#' @description 
#' This function calculates a range of compositional dissimilarity/distance measures 
#' between sites (within studies) based on taxonomic diversity and for the associated
#' HibDatabase.
#' 
#' @details 
#' As is standard for R6 Classes, the HibDatabase is modified in-place.
#' If effort correction (`hib_effort_correction()`) has previously been conducted, 
#' diversity metrics will be estimated using effort corrected measurements.
#' 
#' @param obj A [`HibDatabase`] object.
#' @param metrics Character vector of dissimilarity/distance metrics. 
#' Metrics available via \code{vegan::vegdist}: "manhattan", "euclidean", "canberra", "clark", 
#' "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", 
#' "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", 
#' "hellinger", "aitchison", "robust.aitchison".
#' Other available metrics: "hypervolume", "cosine", "theta", "baselga_bray_balanced", 
#' "baselga_bray_gradient", "baselga_ruzicka_balanced", "baselga_ruzicka_gradient",
#' "horn_dissim", "normalised_euclidean".
#' Metrics can be provided in upper or lower case - internal processing is done 
#' using lower-case - and as as unambiguous abbreviations.
#' Note: "horn" refers to the "Morisita-Horn" index in the vegan package. 
#' "horn_dissim" refers to (1 - Horn overlap), calculated endogenously within 
#' this package.
#' @param binary Logical vector - should metric be calculated using occurence
#' values (`TRUE`) or abundance (`FALSE`, default). Should be a single entry, or
#' of the same length as 'metrics'. Vectors of other lengths will be recycled to 
#' the length of 'metrics'.
#' @param rel_abund Logical vector - determining whether relative abundance 
#' should be used to calculate the dissimilarity metric, only applicable for some metrics, 
#' e.g., euclidean, hypervolume.
#' Note: many metrics include rescaling of abundances as part of the calculation. 
#' Default values is `TRUE`.
#' @param rescale Logical vector - determining if dissimilarity measure should 
#' be rescaled relative to its maximum value - applies to metrics that are not 
#' on a 0-1 scale by default - hypervolume, euclidean (if relative abundance), 
#' chord, hellinger.
#' Default values is `TRUE`.
#' @param baseline_lui Character or NULL, indicating the baseline LUI/s against
#' which compositional dissimilarity is calculated per study. If NULL, compositional 
#' dissimilarity is calculated between all sites within each study.
#' Default values is 'Primary vegetation_Minimal use' - the most pristine LUI within
#' the PREDICTS database.
#' 
#' @return HibDatabase with updated '$diversity_between' tibble
#' 
#' @keywords database, processing
#' 
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
                                 binary = FALSE,
                                 rel_abund = TRUE,
                                 rescale = TRUE,
                                 baseline_lui = "Primary vegetation_Minimal use"
                                 # effort_diff_thr = 0.1
                                 ){

  # TODO:
  # Could transform abund data prior to metric calc - e.g. log1p, log2, sqrt...
  # can affect emphasis of metric - see Fig 5 in below...
  # https://onlinelibrary.wiley.com/doi/full/10.1111/j.1461-0248.2010.01552.x
  
  # Custom metric - 
  # use "custom_metric" in metrics vec
  # add a custom_metric_spec arg - list/named vector...
  # have custom_metric_func, custom_metric_max, custom_metric_name 

  # Assertions
  assertthat::assert_that(
    
    inherits(obj, "HibDatabase"),
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
    
    assertthat::has_name(obj$species, "Taxon_name_entered"),
    assertthat::has_name(obj$species, "Best_guess_binomial"),
    
    assertthat::has_name(obj$sites, "SSBS"),
    assertthat::has_name(obj$sites, "SS"),

    assertthat::has_name(obj$sites, "LUI")
  )
  
  veg_metrics <- c("manhattan", "euclidean", "canberra", "clark", "bray", 
                   "kulczynski", "jaccard", "gower", "altGower", 
                   "morisita", "horn", "mountford", "raup", "binomial", 
                   "chao", "cao", "mahalanobis", "chisq", "chord", 
                   "hellinger", "aitchison", "robust.aitchison")
                   
  # Non-Vegan  
  non_veg_metrics <- c("hypervolume", "cosine", "theta",
                       "baselga_bray_balanced", "baselga_bray_gradient",
                       "baselga_ruzicka_balanced", "baselga_ruzicka_gradient",
                       "horn_dissim", "normalised_euclidean")
  # "baselga_jaccard", "baselga_sorensen")
  
  metrics_ <- c(veg_metrics, non_veg_metrics)
  

  # allow upper/lower-case metrics
  metrics <- tolower(metrics)
  metrics <- unique(metrics)
  # Needs to be altGower
  metrics[metrics == "altgower"] <- "altGower"
  

  # first, combine metrics, binary, rel_abund and rescale into a df
  metrics_spec_df <- data.frame(metric    = metrics,
                                binary    = binary,
                                rel_abund = rel_abund,
                                rescale   = rescale)

  matchable_metrics_idx <- !is.na(pmatch(metrics_spec_df$metric,
                                         metrics_))

  # subset to matchable metrics
  metrics_spec_df <- metrics_spec_df[matchable_metrics_idx,]

  # run match.arg on matchable metrics
  metrics_spec_df$metric <- match.arg(metrics_spec_df$metric, 
                                       metrics_,
                                       several.ok = T)
  
  metrics_spec_df <- unique(metrics_spec_df)

  # %%%
  # if metric is already in hib$diversity_between, ignore
  # TODO: check recalculate...
  if (!is.null(obj$diversity_between)){

    # Need to regen altGower
    metrics_tc <- metrics_spec_df$metric
    
    metrics_tc <- tools::toTitleCase(metrics_tc)
    
    new_metrics_idx <- which(!(metrics_tc %in% obj$diversity_between$Metric))
    # drop from metrics and print the dropped metrics...
    metrics_spec_df <- metrics_spec_df[new_metrics_idx, ]
    dropped_metrics <- metrics_tc[metrics_tc %in% obj$diversity_between$Metric]
    message(crayon::yellow(paste0("\nWarning: the following metrics have already been ",
                                  "added to the HibDatabase and so will not be ",
                                  "recalculated;\n\t",
                                  paste(dropped_metrics, collapse = ", ")
                                  # "\n"
                                  )))
    
    if (nrow(metrics_spec_df) == 0){
        stop("No new biodiversity metrics to calculate.", call. = F)
    }
  }
  # %%%
  
  message(crayon::green("\nPreparing site-contrasts to calculate compositional dissimilarity ...\n"))
  
  
  # Get site comparisons

  # If metrics are binary can also use 'Occurrence' data...
  # If metrics do not use abundance data, these contrasts are filtered out later.
  
  spp_tmp <- obj$species |>
    dplyr::filter(Diversity_metric_type %in% c("Abundance",
                                               "Occurrence")) |>
    
    dplyr::group_by(SS) |>
    dplyr::mutate(# Taxon_name_entered v best_guess_binomial ??
                  n_species = dplyr::n_distinct(Taxon_name_entered)) |>
    dplyr::ungroup() |>
    dplyr::filter(n_species > 1)
  
  if (all(metrics_spec_df$binary == FALSE)){
    spp_tmp <- spp_tmp |>
      dplyr::filter(Diversity_metric_type == "Abundance")

    # print(nrow(spp_tmp))
    if (nrow(spp_tmp) == 0){
      stop("\nError: All metric specifications require Abundance records (binary = FALSE), but none are available.", call. = F)
    }

  }

  
  if (!is.null(baseline_lui)){
    
    # Better appr is to first get site data..., (subset to spp_tmp_SSBS) 
    # then check for baseline LUI...

    # %%%%
    site_tmp <- obj$sites |>
      dplyr::filter(SSBS %in% unique(spp_tmp$SSBS),
                    !is.na(LUI))

    if (nrow(site_tmp) == 0){
      stop("\nError: no site comparisons possible. Check your data and the metric specifications.", call. = F)
    }
    # %%%%%
    # Allow for multiple baseline_lui
    if (!all(baseline_lui %in% unique(site_tmp$LUI))){
      missing_baseline <- baseline_lui[!baseline_lui %in% unique(site_tmp$LUI)]
      baseline_lui <- baseline_lui[baseline_lui %in% unique(site_tmp$LUI)]

      if (length(baseline_lui) > 0){
        message(crayon::yellow(paste0("Warning: specified baseline/s missing, only ",
                                    "using available baseline/s...\n",
                                    "\tMissing baseline/s: ", paste0(missing_baseline, collapse = ", "), "\n",
                                    "\tUsed baseline/s:    ", paste0(baseline_lui, collapse = ", "), "\n"
                                    )))
      } else {
        message(crayon::red(paste0("\nError: no sites avilable with specified baseline LUI/s...\n")))
        stop("Specific baseline not in LUI column.")

      }
      
    }
    assertthat::assert_that(all(baseline_lui %in% unique(site_tmp$LUI)))

    site_tmp$base_lui <- site_tmp$LUI %in% baseline_lui
    # %%%%%
    
    site_tmp <- site_tmp |>
      dplyr::group_by(SS) |>
      dplyr::mutate(n_base_lui = sum(base_lui)) |>
      dplyr::filter(n_base_lui > 0) |>
      dplyr::ungroup() |>
      droplevels()
    # %%%%

    
  } else if (is.null(baseline_lui)){
    site_tmp <- obj$sites |>
      dplyr::filter(SSBS %in% unique(spp_tmp$SSBS),
                    !is.na(LUI)) |>
      droplevels()
  }
  

  n_spp <- spp_tmp |>
    dplyr::select(Best_guess_binomial) |>
    dplyr::n_distinct()
  n_spp_abund <- spp_tmp |>
    dplyr::filter(Diversity_metric_type == "Abundance") |>
    dplyr::select(Best_guess_binomial) |>
    dplyr::n_distinct()

  
  # TODO: set any NA measures to 0 for safety
  if (sum(is.na(spp_tmp$Measurement)) > 0){
    message(crayon::yellow(paste0("\nWarning: NAs detected in 'Measurement'.  Setting NAs to zeros.")))
    spp_tmp$Measurement[is.na(spp_tmp$Measurement)] <- 0
  }

  assertthat::assert_that(
    sum(is.na(spp_tmp$Measurement)) == 0,
    is.numeric(spp_tmp$Measurement)
  )
  
  # Need to get combinations of data... - s1, s2, ...
  # %%%
  # From https://adrianadepalma.github.io/BII_tutorial/bii_example.html
  
  # get a vector of each study to loop over
  studies <- unique(as.character(site_tmp$SS))
  
  site_comp_df <- purrr::map_dfr(.x = studies, .f = function(x){
    
    # pull out the sites that are the baseline lui (we only want to use comparisons with the baseline)
    if (!is.null(baseline_lui)){
      site_data <- dplyr::filter(site_tmp, SS == x) |>
        # pull out the SSBS and LandUse information
        dplyr::select(SSBS, base_lui)
      
      baseline_sites <- site_data |>
        dplyr::filter(base_lui == 1) |>
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
      dplyr::rename(s1 = Var1, s2 = Var2) |>
      # make the values characters rather than factors
      dplyr::mutate(s1 = as.character(s1),
                    s2 = as.character(s2)) |>
      # remove the comparisons where the same site is being compared to itself
      dplyr::filter(s1 != s2) |>      
                    # add the full name
      dplyr::mutate(contrast = paste(s1, s2, sep = "_"),
                    # add the study id
                    SS = as.character(x))
    
    return(site_comparisons)
  })
  
  # Add effort for ssbs1/2 ??

  # if nrow(site_comp_df == 0) stop...
  if (nrow(site_comp_df) == 0){
    message(crayon::red(paste0("\nError: no site comparisons possible...\n\n",
                               "Maybe check that the specified baseline LUI/s ('", 
                               paste(baseline_lui, collapse = ", "),
                               "') is present ...\n"))
    )
    stop("Unable to generate any site comparisons, see above for details.", call. = FALSE)
  }


  message(crayon::green(paste0(nrow(site_comp_df), 
                               " site comparisons to make ..."))
          )
  

  # add Diversity_metric_type to site_comp_df
  site_comp_df <- site_comp_df |>
    dplyr::left_join(site_tmp[, c("SSBS", "Diversity_metric_type")] |>
                        dplyr::rename(Diversity_metric_type_s1 = Diversity_metric_type),
                      by = c("s1" = "SSBS")) |>
    dplyr::left_join(site_tmp[, c("SSBS", "Diversity_metric_type")] |>
                        dplyr::rename(Diversity_metric_type_s2 = Diversity_metric_type),
                      by = c("s2" = "SSBS")) 


  # assert that they are same for s1/s2
  assertthat::assert_that(
    all(site_comp_df$Diversity_metric_type_s1 == site_comp_df$Diversity_metric_type_s2)
  )

  site_comp_df <- site_comp_df |>
    dplyr::select(-Diversity_metric_type_s2) |>
    dplyr::rename(Diversity_metric_type = Diversity_metric_type_s1)



  # Some checks of metric specs and diversity data types available
  
  un_div <- unique(site_comp_df$Diversity_metric_type) |> as.character()
  # if all binary are false, / 
  # if some binary are false, but only occurrence in site_comp_df, warning and drop metrics
  # if all binary and only only occurrence in site_comp_df - drop all metrics warning/stop
  if ((all(metrics_spec_df$binary == F)) &
      (!any("Abundance" %in% un_div)) ){

    stop("\nError: All metric specifications require Abundance records (binary = FALSE), but none are available.", call. = F)
  }

  if ((any(metrics_spec_df$binary == F)) &
      (!any("Abundance" %in% un_div)) ){

    old_metr <- metrics_spec_df$metric
    metrics_spec_df <- metrics_spec_df |>
      dplyr::filter(binary)

    new_metr <- metrics_spec_df$metric
    dropped_metr <- setdiff(old_metr, new_metr)

    message(crayon::yellow(paste0("\nWarning: the following metrics have been dropped", 
                                  " as the specifications require Abundance records",
                                  " (binary = FALSE), but none are available; \n\t",
                                  paste(tools::toTitleCase(dropped_metr), 
                                       collapse = ", ")  )))
    rm(old_metr, new_metr, dropped_metr)
  }


  # split to veg and non_veg spec dfs
  veg_metrics_df <- subset(metrics_spec_df,
                           metric %in% veg_metrics)
  non_veg_metrics_df <- subset(metrics_spec_df,
                               metric %in% non_veg_metrics)

  message(crayon::green("\nCalculating compositional dissimilarity ...\n"))
  

  # do vegan-based metrics...
  # %%%
  veg_df <- NULL
  
  if (nrow(veg_metrics_df) > 0){

    message(crayon::white(paste0("\tVegan-based dissimilarity metrics (",
                                 paste(tools::toTitleCase(veg_metrics_df$metric), 
                                       collapse = ", "),
                                 ") ...\n")))


    veg_df <- purrr::map_dfr(.x = studies,
                                 .f = function(ss = .x, 
                                               spp_data = spp_tmp,
                                               comp_df = site_comp_df,
                                               metric_df = veg_metrics_df){
                                   
                                   spp_site_data <- spp_data |>
                                     # filter out the SS that matches the SS we're interested in
                                     dplyr::filter(SS == ss) |>
                                     # pull out the site name, species name and abundance information
                                     # pass to hib_calc_veg...
                                     dplyr::select(SSBS, Taxon_name_entered, Measurement) |>                                  
                                     # pivot so that each column is a species and each row is a site
                                     tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement,
                                                        values_fill = 0) |>
                                     
                                     # set the rownames to the SSBS and then remove that column
                                     tibble::column_to_rownames("SSBS")
                                   
                                   # Get Diversity metric_type of each SSBS
                                   
                                   div_metr_type <- comp_df |>
                                     dplyr::filter(SS == ss) |>
                                     dplyr::pull(Diversity_metric_type) |>
                                     unique() |> as.character()

                                   
                                   assertthat::assert_that(
                                     length(div_metr_type) == 1
                                   )

                                   # TODO: check at least 2 species are present (i.e 2 cols or more)
                                   # Depending on metric_vec values, calc metrics of interest
                                   # %%%
                                   out_df <- purrr::pmap(.l = metric_df,
                                                            .f = function(metric,
                                                                          binary,
                                                                          rel_abund,
                                                                          rescale,
                                                                          spp_data = spp_site_data,
                                                                          div_metric_type = div_metr_type){
                                                              hib_calc_veg(spp_data, 
                                                                           metric,
                                                                           binary,
                                                                           rel_abund,
                                                                           rescale,
                                                                           div_metric_type)
                                                            }) |>
                                     dplyr::bind_rows()
                                   # %%%
                                   
                                   out_df <- out_df |>
                                     dplyr::filter(contrast %in% comp_df$contrast)
                                   
                                   out_df <- tidyr::pivot_wider(out_df,
                                                                names_from = metric_nm, 
                                                                values_from = Distance)
                                   out_df$SS <- ss
                                   
                                   return(out_df)
                                 })
    
    veg_df <- veg_df |>
      dplyr::select(-contrast)
  }
  # %%%
  
  # do other, non-vegan metrics...
  non_veg_df <- NULL
  
  if (nrow(non_veg_metrics_df) > 0){
    
    message(crayon::white(paste0("\tNon-vegan dissimilarity metrics (",
                                 paste(tools::toTitleCase(non_veg_metrics_df$metric), 
                                       collapse = ", "),
                                 ") ...\n")))
    
    non_veg_df <- purrr::map2_dfr(.x = site_comp_df$s1,
                              .y = site_comp_df$s2,
                              .f = function(s1 = .x, 
                                            s2 = .y, 
                                            spp_data = spp_tmp,
                                            metric_d = non_veg_metrics_df){

                                spp_site_data <- spp_data |>
                                  # filter out the SSBS that matches the pair of sites we're interested in
                                  dplyr::filter(SSBS %in% c(s1, s2)) |>
                                  # pull out the site name, species name and abundance information
                                  dplyr::select(SSBS, Taxon_name_entered, Measurement) |>
                                  # pivot so that each column is a species and each row is a site
                                  tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement,
                                                     values_fill = 0) |>
                                  # set the rownames to the SSBS and then remove that column
                                  tibble::column_to_rownames("SSBS")
                                # TODO: Set row 1 to s1 (i.e. baseline, or control site)
                                # Only needed if metric is directional...

                                ss <- spp_data |>
                                  # filter out the SSBS that matches the pair of sites we're interested in
                                  dplyr::filter(SSBS %in% c(s1, s2)) |>
                                  dplyr::pull(SS) |>
                                  unique()
                                
                                out_df <- data.frame(SSBS1 = s1, 
                                                     SSBS2 = s2,
                                                     SS    = ss)

                                # Get Diversity metric_type
                                div_metr_type <- spp_data |>
                                     dplyr::filter(SS == ss) |>
                                     dplyr::pull(Diversity_metric_type) |>
                                     unique() |> as.character()

                                assertthat::assert_that(
                                  length(div_metr_type) == 1
                                )
                                
                                # Depending on sim_metric values, calc metrics of comp-sim
                                # non_veg_metrics_df
                                out_met <- purrr::pmap(.l = non_veg_metrics_df,
                                                          .f = function(metric,
                                                                        binary,
                                                                        rel_abund,
                                                                        rescale,
                                                                        spp_data = spp_site_data,
                                                                        div_metric_type = div_metr_type){
                                                            hib_calc_nonveg(spp_data, 
                                                                            metric,
                                                                            binary,
                                                                            rel_abund,
                                                                            rescale,
                                                                            div_metric_type)
                                                          }) |>
                                  dplyr::bind_cols()
                                
                                out_df <- out_df |> dplyr::bind_cols(out_met)
                                return(out_df)
                              }
    )
  }
  
  # combine veg and non_veg
  # s1, s2, ...
  if ((!is.null(veg_df)) & (!is.null(non_veg_df))){
    veg_df <- dplyr::left_join(veg_df, non_veg_df,
                                by = c("SSBS1", "SSBS2", "SS")) 
  }
  if (is.null(veg_df) & (!is.null(non_veg_df))){
    veg_df <- non_veg_df
  }
  
  # add processing flags per metric -  

  # Calc is_empty flags once and add to each sub df?
  empty_sites1 <- spp_tmp |>
    dplyr::group_by(SSBS) |>
    dplyr::select(SSBS, Taxon_name_entered, Measurement) |>
    # pivot so that each column is a species and each row is a site
    tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement,
                       values_fill = 0) |>
    tibble::column_to_rownames("SSBS") 
    # %%%%%

  empty_sites1$SSBS1_is_empty <- rowSums(empty_sites1, na.rm = T) == 0

  empty_sites1 <- empty_sites1 |>
    tibble::rownames_to_column("SSBS1") |>
    dplyr::select(SSBS1, SSBS1_is_empty)
    # %%%%%
    
  
  empty_sites2 <- empty_sites1 |>
    dplyr::rename(SSBS2 = SSBS1,
                  SSBS2_is_empty = SSBS1_is_empty)
  
    
  veg_df <- veg_df |>
    dplyr::left_join(empty_sites1, by = "SSBS1") |>
    dplyr::left_join(empty_sites2, by = "SSBS2")


  veg_df <- veg_df |>
    dplyr::left_join(site_comp_df[, c("s1", "s2", "Diversity_metric_type")] |>
                        dplyr::rename(SSBS1 = s1,
                                      SSBS2 = s2), 
                    by = c("SSBS1", "SSBS2"))
  # %%%%%
  # 
  metrics_out <- colnames(veg_df)
  metrics_out <- metrics_out[!metrics_out %in% c("SSBS1", "SSBS2", "SS",
                                                 "Diversity_metric_type",
                                                 "SSBS1_is_empty", "SSBS2_is_empty")] 

  metrics_out_ls <- purrr::map(metrics_out,
                               function(x,
                                        df = veg_df){
                                 return(df[,c("SSBS1", "SSBS2", "SS", x,
                                              "SSBS1_is_empty", "SSBS2_is_empty",
                                              "Diversity_metric_type")])
                               })
  
  out_tib <- dplyr::tibble(Metric = metrics_out,
                           Data = metrics_out_ls,
                           N_species_total = n_spp
                           )
  # Add Flag based on metric and specific flags...
  # need to match sentence-case metric names to lower-case metric names...
  metrics_spec_df <- metrics_spec_df |>
    dplyr::mutate(Metric = tools::toTitleCase(metric))
  
  metrics_spec_df <- metrics_spec_df |>
    dplyr::mutate(baseline_lui = paste(baseline_lui, collapse = ", "))
  if (is.null(baseline_lui)){
    metrics_spec_df <- metrics_spec_df |>
      dplyr::mutate(baseline_lui = NA)
  }
  
  # Improve column names...
  metrics_spec_df_ <- metrics_spec_df[,c("Metric", "binary", "rel_abund", "rescale",
                                         "baseline_lui")]
  colnames(metrics_spec_df_) <- c("Metric", "Binary", "Relative_abundance", "Rescale",
                                  "Baseline_LUI")
  out_tib <- out_tib |>
    dplyr::left_join(metrics_spec_df_, by = "Metric")
  

  # %%%%%
  # If Binary is False, filter internal dfs to Abundance only..
  for (i in 1:nrow(out_tib)){
    if (out_tib$Binary[i] == FALSE){
      out_tib$Data[[i]] <- out_tib$Data[[i]] |>
                            dplyr::filter(Diversity_metric_type == "Abundance")

      # For non-binary distances, set n_spp to n_spp_abund
      out_tib$N_species_total[i] <- n_spp_abund
    }
  }
  

  # %%%%%
  # add in max_distance/max_values
  # Based on veg and non-veg...

  # veg
  # if(method %in% c(3,4,5,7,8,10,11,13,17)) # index in 0..1
  #           1
  #       else if (method %in% c(19,20)) # chord, hellinger
  #           sqrt(2)
  #       else if (method == 9) # Mountford
  #           log(2)
  #       else # no fixed upper limit
  #           NA
  
  veg_max_df <- data.frame(Metric = tools::toTitleCase(c("manhattan", "euclidean", "canberra", "bray", # 4
                                                        "kulczynski", "gower", "morisita", "horn", #8
                                                        "mountford", "jaccard", "raup", "binomial", "chao", #13
                                                        "altGower", "cao", "mahalanobis", "clark", "chisq", "chord", #19
                                                        "hellinger", "aitchison", "robust.aitchison")),
                          Theoretical_max = c(NA, NA, 1, 1, 
                                             1,  1, 1, 1,
                                             log(2), 1, 1, NA, 1,
                                             NA, NA, NA, 1, NA, sqrt(2), 
                                             sqrt(2), NA, NA))
  # only if Euclidean is rel_abund
  if ("Euclidean" %in% metrics_spec_df_$Metric){
    euc_rel <- metrics_spec_df_$Relative_abundance[metrics_spec_df_$Metric == "Euclidean"]
    if (euc_rel == T){
      veg_max_df$Theoretical_max[veg_max_df$Metric == "Euclidean"] <- sqrt(2)
    }
  }

  # Max gower dist for 2 sites is 1, rel_abund or not...
  # ...  sum (diff_i / diff_i) / N
  # if diff_i != 0, ratio is 1, 
  # returns prop of spp with different abundances
  # 

  # non-veg
  # c("hypervolume", "cosine", "theta",
                        # "baselga_bray_balanced", "baselga_bray_gradient", 
                        # "baselga_ruzicka_balanced", "baselga_ruzicka_gradient")
  # hypervol = NA, or 1 if rel_abund
  # cosine = 0-1 (set 1 as max, convert to dissim)
  # theta  = pi/2
  # all others, 0-1 - check with baselga...
  # normalised_euclidean - 1.

  non_veg_max_df <- data.frame(Metric = tools::toTitleCase(c("hypervolume", "cosine", "theta",
                                                             "baselga_bray_balanced", "baselga_bray_gradient", 
                                                             "baselga_ruzicka_balanced", "baselga_ruzicka_gradient",
                                                             "horn_dissim", "normalised_euclidean")),
                              Theoretical_max = c(NA, 1, pi/2,
                                                 1, 1,
                                                 1, 1,
                                                 1, 1))
  # only if Hypervolume is rel_abund
  if ("Hypervolume" %in% metrics_spec_df_$Metric){
    vol_rel <- metrics_spec_df_$Relative_abundance[metrics_spec_df_$Metric == "Hypervolume"]
    if (vol_rel == T){
      non_veg_max_df$Theoretical_max[non_veg_max_df$Metric == "Hypervolume"] <- 1
    }
  }


  max_df <- dplyr::bind_rows(veg_max_df,
                             non_veg_max_df)

  out_tib <- out_tib |>
    dplyr::left_join(max_df, by = "Metric") |>
    dplyr::mutate(Rescale = dplyr::case_when(is.na(Theoretical_max) ~ FALSE,
                                              T ~ TRUE))

  # %%%%%
  
  div_empty <- is.null(obj$diversity_between)
  
  if (!div_empty){
    obj$diversity_between <- dplyr::bind_rows(obj$diversity_between, out_tib)
  }  
  if (div_empty){
    obj$diversity_between <- out_tib
  }
  
  obj$last_update <- lubridate::now()

  return(obj)
}

