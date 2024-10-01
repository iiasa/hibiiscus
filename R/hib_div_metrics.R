# TODO: recalculate arg - recalculates metrics specified - over-rules dropping of
# previously calculated metrics - e.g if re-filtering DB...
# TODO: check whether or not to use N_spp calculations per site/total 
# Best_guess_binomial isn't always useful in PREDICTS - e.g. if only a higher 
# taxonomic level is recorded...
# TODO: add unequalness 
# TODO: Hill numbers of various powers? -  
# TODO: check limits of metric outputs: none should be negative, many should be <=1, 
# others have a specific max

#' Calculate number of species (best guess) and original taxa names
#' 
#' @description 
#' This function calculates number of species contributing to diversity metrics.
#' @details 
#' Expects a subset ofthe species data.frame from a [`HibDatabase`] object.
#' @param df A subset of the species data.frame from a [`HibDatabase`] object.
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   out <- get_n_spp_taxa(df)
#' } 
# @export
#' @noRd
get_n_spp_taxa <- function(df, measure_var){
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(df, "Taxon_name_entered"),
    assertthat::has_name(df, "Best_guess_binomial"),
    assertthat::has_name(df, measure_var)
  )

  df <- df |>
    dplyr::filter(!!dplyr::sym(measure_var) > 0)

  n_spp <- df |>
    dplyr::filter(Best_guess_binomial != '' &
                  !is.na(Best_guess_binomial)) |>
    dplyr::select(Best_guess_binomial) |>
    dplyr::n_distinct()

  n_orig_taxa <- df |>
    dplyr::filter(Taxon_name_entered != '' &
                  !is.na(Taxon_name_entered)) |>
    dplyr::select(Taxon_name_entered) |>
    dplyr::n_distinct()

  return(list(N_spp = n_spp,
              N_taxa = n_orig_taxa))
}


#' Calculate total abundance
#' 
#' @description 
#' This function calculates total abundance for each site measuring abundance.
#' @details 
#' Expects the species data.frame from a [`HibDatabase`] object.
#' Total_abundance = sum(x[i]), where x[i] is the abundance of species i within 
#' a given site.
#' Rescaled_total_abundance = Total_abundance[j,k]/max(Total_abundance[j,k]), where
#' the total abundance in eack site j is divided by the maximum total abundance 
#' in study k.
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within hib_calc_div_within().
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_ta(df, "Measurement")
#' } 
# @export
#' @noRd
hib_calc_ta <- function(spp_df, measure_var){
  
  message(crayon::white("\tTotal abundance ...\n"))
  
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), # obj$species
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
                  # Rescaled_total_abundance = Total_abundance / max(Total_abundance)) |>
    dplyr::ungroup() |>
    dplyr::select(-max_ta)
  
  # Ensure all rescaled values between 0 and 1! - avoid rounding errors
  abund_df$Rescaled_total_abundance[abund_df$Rescaled_total_abundance < 0] <- 0
  abund_df$Rescaled_total_abundance[abund_df$Rescaled_total_abundance > 1] <- 1
  
  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  # n_guess_taxa <- sub_df |>
  #   dplyr::select(Best_guess_binomial, Taxon_name_entered) |>
  #   dplyr::mutate(Best_guess_taxon = Best_guess_binomial) |>
  #   dplyr::mutate(Best_guess_taxon = dplyr::case_when(is.na(Best_guess_taxon) ~ Taxon_name_entered,
  #                                                           Best_guess_taxon == ""   ~ Taxon_name_entered,
  #                                                           T ~ Best_guess_taxon))|>
  #   dplyr::select(Best_guess_taxon) |>
  #   dplyr::n_distinct()

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
#' @details 
#' Expects the species data.frame from a [`HibDatabase`] object.
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within hib_calc_div_within().
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_sr(spp_df, measure_var)
#' } 
# @export
#' @noRd
hib_calc_sr <- function(spp_df, measure_var){
  message(crayon::white("\tSpecies richness ...\n"))
  
  # Assertion checks
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


#' Calculate simpson's diversity
#' 
#' @description 
#' This function calculates simpson's diversity per site.
#' @details 
#' Expects the species data.frame from a [`HibDatabase`] object.
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within hib_calc_div_within().
#' @param msg A [`logical`], indicating if default progress message should be 
#' displayed or not.
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_simpson(df, "Measurement")
#' } 
#' @noRd
hib_calc_simpson <- function(spp_df, measure_var, msg = TRUE){
  # optional message - drop if used for Hill-Simpson
  assertthat::assert_that(
    is.logical(msg)
  )
  if (msg){
    message(crayon::white("\tSimpson diversity ...\n")) # Calculating within-site   
  }
  
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), # obj$species
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  # Filter to useable records
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance") 

  simp_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::summarise(Simpson_diversity = 1 - sum( ( !!dplyr::sym(measure_var) / sum(!!dplyr::sym(measure_var)) )^2 ) ) |>
    dplyr::ungroup()
  # %%%%%
  # Ensure all values between 0 and 1!
  simp_df$Simpson_diversity[simp_df$Simpson_diversity < 0] <- 0
  simp_df$Simpson_diversity[simp_df$Simpson_diversity > 1] <- 1
  
  # TODO: above is gini-simpson?/uncorrected PIE....
  # Above is 'Simpson's Index of Diversity' - 0 = no div, 1 = max div

  # Std Simpson index is...
  # sum( ( !!dplyr::sym(measure_var) / sum(!!dplyr::sym(measure_var)) )^2 )
  # 1 = no div, 0 = max div

  # KEY Q
  # TODO: what if total abundance == 0 - Simpson will be undefined - should it be set to 0??
  # But, if only 1 species, index value is 0 - i.e no diversity if all one species...

  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  out_tib <- dplyr::tibble(Metric = "Simpson_diversity",
                           Data = list(simp_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa) 

  return(out_tib)
}


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
#' Expects the species data.frame from a [`HibDatabase`] object.
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within hib_calc_div_within().
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   tib <- hib_calc_chao1(spp_df, measure_var)
#' } 
#' @noRd
hib_calc_chao1 <- function(spp_df, measure_var){
  message(crayon::white("\tChao1 species richness ...\n"))
  
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), # obj$species
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Diversity_metric_is_suitable_for_Chao"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance" &
                    Diversity_metric_is_suitable_for_Chao) |>
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
#' Expects the species data.frame from a [`HibDatabase`] object.
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within hib_calc_div_within().
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   out <- hib_calc_rsr(spp_df, measure_var)
#' } 
#' @noRd
hib_calc_rsr <- function(spp_df, measure_var){
  message(crayon::white("\tRarefied species richness ...\n"))
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), # obj$species
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Diversity_metric_is_suitable_for_Chao"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  rarefy_comm <- function(abund_distr, N){
    comm <- rep(letters[1:length(abund_distr)], abund_distr)
    samp <- sample(comm, N)
    sr <- length(table(samp))
    # <- length(unique(samp)) # more efficient?
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
  
  pb <- utils::txtProgressBar(min = 1, max = 100, style = 3)
  
  rsr_df <- purrr::map_dfr(seq(1,100),
                           function(x, df = rsr_df){
                             tmp <- df |>
                               dplyr::group_by(SSBS, SS) |>
                               dplyr::summarise(sr = rarefy_comm(!!dplyr::sym(measure_var), 
                                                                 min_abund[1]),
                                                min_abund = min_abund[1],
                                                iter = x)
                             
                             utils::setTxtProgressBar(pb, x)
                             # print(head(tmp))
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



# %%%%%
#' Calculate arithmetic mean abundance
#' 
#' @description 
#' This function calculates arithmetic mean abundance for each site measuring abundance.
#' @details 
#' Expects the species data.frame from a [`HibDatabase`] object.
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within hib_calc_div_within().
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   out <- hib_calc_ma(spp_df, mv)
#' } 
#' @noRd
hib_calc_ma <- function(spp_df, measure_var){
  message(crayon::white("\tMean abundance ...\n"))
  
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), # obj$species
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  

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

    # dplyr::mutate(Rescaled_mean_abundance = Mean_abundance / max(Mean_abundance)) |>
    # dplyr::ungroup()
  
  # Ensure all rescaled values between 0 and 1!
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


# %%%%%
#' Calculate geometric mean abundance
#' 
#' @description 
#' This function calculates geometric mean abundance for each site measuring abundance.
#' @details 
#' Expects the species data.frame from a [`HibDatabase`] object.
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within hib_calc_div_within().
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   out <- hib_calc_gma(spp_df, measure_var)
#' } 
#' @noRd
hib_calc_gma <- function(spp_df, measure_var){
  message(crayon::white("\tGeometric mean abundance ...\n"))
  
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), # obj$species
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance")

  abund_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::filter(!!dplyr::sym(measure_var) > 0) |>
    # dplyr::summarise(Geometric_mean_abundance = dplyr::case_when(dplyr::n() == 0 ~ 0,
    #                                                              T ~ exp(mean(log(!!dplyr::sym(measure_var)))) ) ) |>
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
    # if max Mean_abundance is 0..., leave as is..
    dplyr::mutate(max_gma = max(Geometric_mean_abundance),
                  Rescaled_geometric_mean_abundance = dplyr::case_when(max_gma == 0 ~ Geometric_mean_abundance,
                                                                       T ~ Geometric_mean_abundance / max_gma)) |>
    dplyr::ungroup() |>
    dplyr::select(-max_gma)
    # dplyr::mutate(Rescaled_geometric_mean_abundance = Geometric_mean_abundance / max(Geometric_mean_abundance)) |>
    # dplyr::ungroup()
  
  # TODO: ensure all rescaled values between 0 and 1!
  # abund_df$Rescaled_geometric_mean_abundance[abund_df$Rescaled_geometric_mean_abundance < 0] <- 0
  # abund_df$Rescaled_geometric_mean_abundance[abund_df$Rescaled_geometric_mean_abundance > 1] <- 1
  
  # Get N species used
  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  out_tib <- dplyr::tibble(Metric = "Geometric_mean_abundance",
                           Data = list(abund_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa)
  
  return(out_tib)
}


# %%%%%
#' Calculate shannon diversity
#' 
#' @description 
#' This function calculates shannon diversity for each site measuring abundance
#' @details 
#' Expects the species data.frame from a [`HibDatabase`] object.
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within hib_calc_div_within().
#' @param msg A [`logical`], indicating if default progress message should be 
#' displayed or not.
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   out <- hib_calc_shannon(spp_df, mv)
#' } 
#' @noRd
hib_calc_shannon <- function(spp_df, measure_var, msg = TRUE){
  assertthat::assert_that(
    is.logical(msg)
  )
  if (msg){
    message(crayon::white("\tShannon diversity ...\n"))  
  }
  
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), # obj$species
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance")

  div_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::mutate(p_i = !!dplyr::sym(measure_var) / sum(!!dplyr::sym(measure_var)) ) |>
    # If all p_i == 0, or NA, Shannon_diversity = NA
    dplyr::summarise(Shannon_diversity =  dplyr::case_when(all(p_i == 0) ~ NA,
                                                          all(is.na(p_i)) ~ NA,
                                                           T ~ -1 * sum( p_i * log(p_i), na.rm = TRUE ) ),
    # dplyr::summarise(Shannon_diversity =  -1 * sum( p_i * log(p_i), na.rm = TRUE ), #) |> 
                    logN = log(sum(!!dplyr::sym(measure_var) > 0)) ) |>
    dplyr::ungroup()

  # Get N species used
  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  # Store in obj$diversity
  out_tib <- dplyr::tibble(Metric = "Shannon_diversity",
                           Data = list(div_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa)
  
  return(out_tib)
}



# %%%%%
# TODO: Hulbert's PIE
# Hurlbert's (1971) Probability of Interspecific Encounter (PIE):
# PIE = (N/(N-1)) * (1 - sum(p^2))
# where:
# p = the proportion a specific taxon composes of a sample
# N = the total abundance of organisms
#
# This measurement is equivalent to Simpson's Index but includes a correction
# factor base on the total abundance of organisms.

#' Calculate PIE (Hulbert's (1971) Probability of Interspecific Encounter)
#' 
#' @description 
#' This function calculates PIE per site.
#' @details 
#' Expects the species data.frame from a [`HibDatabase`] object.
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within hib_calc_div_within().
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   out <- hib_calc_pie(df, mv)
#' } 
#' @noRd
hib_calc_pie <- function(spp_df, measure_var){
  message(crayon::white("\tProbability of Interspecific Encounter (PIE) ...\n"))
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), # obj$species
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
  spp_df <- spp_df |>
    dplyr::filter(Diversity_metric_type == "Abundance") 

  pie_df <- spp_df |>
    dplyr::group_by(SSBS, SS) |>
    dplyr::summarise(adj_val = sum( !!dplyr::sym(measure_var) ) / ( sum( !!dplyr::sym(measure_var) ) - 1 ),
                     PIE_ = 1 - sum( ( !!dplyr::sym(measure_var) / sum(!!dplyr::sym(measure_var)) )^2 ),
                     PIE = adj_val * PIE_,
                     sPIE = (1 / (1 - PIE)) ) |>
    dplyr::ungroup() |>
    dplyr::select(SSBS, SS, PIE, sPIE)
  
  # Ensure all PIE values between 0 and 1!
  pie_df$PIE[pie_df$PIE < 0] <- 0
  pie_df$PIE[pie_df$PIE > 1] <- 1
  
  # Get N species used
  spp_info <- get_n_spp_taxa(spp_df, measure_var)

  out_tib <- dplyr::tibble(Metric = "PIE",
                           Data = list(pie_df),
                           N_species_total = spp_info$N_spp,
                           N_taxa_orig_total = spp_info$N_taxa)
  
  return(out_tib)
}



# %%%%%
#’ On the diversity-based measures of equalness and evenness
#’
#’ @param f numeric vector with absolute frequencies
#’
#’ @param q order of Hill numbers
#’
#’ @export
#’ @references "H-R Gregorius, EM Gillet. ####.
#’ On the diversity-based measures of equalness and evenness.
#’ Methods in Ecology and Evolution,
#’ https://doi.org/####.###"
#’
#’ @returns
#’
#’ @examples
#’ f = c( 10, 20, 30, 40, 50, 60 )
#’ q = 1
#’ e( f, q )
get_eq_ev = function(f, q = 1){
  
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

  # out <- data.frame(N = n,
  #                   div_freq = vf,
  #                   div_step = vh,
  #                   ev = ef,
  #                   eq = eh)
  out <- eh
  return(out)
}


#' Calculate equalness 
#' 
#' @description 
#' This function calculates 'equalness' per site, following Gregorius and Gillet, 
#' 2024.
#' @details 
#' Expects the species data.frame from a [`HibDatabase`] object.
#' @param spp_df The species data.frame from a [`HibDatabase`] object.
#' @param measure_var Character, either 'Measurement' or 'Effort_corrected_measurement'
#' generated within hib_calc_div_within().
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_equalness(df, "Measurement")
#' } 
#' @noRd
hib_calc_equalness <- function(spp_df, measure_var){
  message(crayon::white("\tEqualness ...\n")) # Calculating within-site   
  
  # NOTE: could use evenness from Gregorious? - similar to mine, but - 1 in 
  # numerator and denominator
  
  # Assertion checks
  assertthat::assert_that(
    assertthat::has_name(spp_df, "SS"), # obj$species
    assertthat::has_name(spp_df, "SSBS"),
    assertthat::has_name(spp_df, measure_var),
    assertthat::has_name(spp_df, "Diversity_metric_type"),
    assertthat::has_name(spp_df, "Taxon_name_entered"),
    assertthat::has_name(spp_df, "Best_guess_binomial")
  )
  
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
# %%%%%


#' Calculate within-site diversity metrics
#' 
#' @description 
#' This function calculates diversity per site, wrapper for the above functions
#' @details 
#' Expects a [`HibDatabase`] object.
#' @param obj A [`HibDatabase`] object.
#' @param metrics A [`character`] vector, any/all of `"Total_abundance"`, "Species_richness", 
#' "Chao1_diversity",  "Simpson_diversity", "Rarefied_species_richness",
#' "Mean_abundance", "Geometric_mean_abundance", "Shannon_diversity",
#' "HillShannon_diversity", "HillSimpson_diversity", "PIE". 
#' Metrics can be provided in upper or lower case - internal processing is done 
#' using lower-case - and as unambiguous abbreviations.
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   obj <- hib_calc_div_within(obj, c("Total_abundance", "Species_richness"))
#' } 
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
  # message(paste(metrics, collapse = ", "))


  # %%%
  # TODO: if metric is already in hib$diversity_within,
  if (!is.null(obj$diversity_within)){
    metrics_tc <- tools::toTitleCase(metrics)
    
    # %%%
    # TODO - "HillShannon_diversity", "HillSimpson_diversity", "PIE" need to be
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
    if (lenght(dropped_metrics) > 0){
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
                           # tmp <- hib_calc_shannon(db_obj, msg = FALSE)
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
                         # Effective N spp (Shannon) - exp(Shannon) ?? - Hill-Shannon
                         # Effective N spp (Simpson) - 1/Simpson
                         # Coverage-based richness - N spp if 80% of comm sampled
                         
                         # %%%
                         # Hill-Shannon: exp(shannon)
                         if (metric == "hillshannon_diversity"){
                           message(crayon::white("\tHill-Shannon diversity ...\n"))
                           # tmp <- hib_calc_shannon(db_obj, msg = FALSE)
                           tmp <- hib_calc_shannon(spp_df, mv, msg = FALSE)
                           tmp$Metric <- "HillShannon_diversity"
                           tmp$Data[[1]]$HillShannon_diversity <- exp(tmp$Data[[1]]$Shannon_diversity)
                           
                           tmp$Data[[1]] <- tmp$Data[[1]] |>
                             dplyr::select(SSBS, SS, HillShannon_diversity)
                           return(tmp)
                         }
                         
                         # Hill-Simpson: 1/simpson 
                         # TODO: why is my Simpson, 1 - x... - i.e. PIE
                         if (metric == "hillsimpson_diversity"){
                           message(crayon::white("\tHill-Simpson diversity ...\n"))
                           # tmp <- hib_calc_simpson(db_obj, msg = FALSE)
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
# %%%%%



# vegan R package can calc many "distance" metrics in an approach analagous
# to the dist() function, may be much quicker, includes Horn-Morista...
# Takes a site (rows) by spp (columns) matrix as input

# vegan metrics:
# c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", 
# "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", 
# "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", "robust.aitchison")

# conflicts between current metrics and vegan metrics...
# c("jaccard", "bray")

# what about sorensen - what is this?? - binary bray...
# Sorensen similarity presented by Santini et al 2018 is equivalent to 1 - bray,
# when bray is calculated by vegan::vegdist()

# ruzika - equivalent to quantitative jaccard??

# also need to record if one site has 0 total abundance...

# use_relative_abundance should only apply to euclidean and hypervolume.. - maybe bray..
# all other metrics specifically include some sort of scaling


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
#' @param div_metric_type A [`character`] describing the type of biodiversity records
#' at the site of interest - one of "Occurrence" or "Abundance". Set by hib_cacl_div_between.
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
                         rescale,
                         div_metric_type
                         ){
  
  assertthat::assert_that(
    is.character(metric),
    is.logical(binary) || is.null(binary),
    is.logical(rel_abund) || is.null(rel_abund),
    is.logical(rescale) || is.null(rescale),
    is.character(div_metric_type)
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
  
  # message(metric)
  # if ((rel_abund) & (metric == "euclidean")){
  #   row_abund <- rowSums(spp_site_data)
  #   row_abund[row_abund==0] <- 1
  #   spp_site_data <- spp_site_data/row_abund
  # }
  # TODO: NA checks...
  # if binary
  # # if !binary, filter to Abundance Diversity_metric_type - then drop col

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
  
  # TODO: if binary is false, and div_metr_type is Occurrence
  # Set to NAs?

  return(dist_df)
}



#' Calculate cosine similarity metric
#' 
#' @description 
#' This function calculates cosine similarity between the relative abundances of 
#' species at two sites
#' @details 
#' Expects two vectors, one each for the abundance/occurence of species at two sites
#' assuming that values at the same index (e.g., x1[i], x2[i]) correspond to the 
#' same species/taxonomic entity
#' @param x1 A [`numeric`] vector for abundance/occurence of species at site 1
#' @param x2 A [`numeric`] vector for abundance/occurence of species at site 2
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   cosine_sim <- calc_cosine_sim(c(0,1,1), c(1,0,1))
#' } 
#' @noRd
calc_cosine_sim <- function(x1, x2){
  assertthat::assert_that(
    length(x1) == length(x2)
  )

  cosine <- sum( x1 * x2 , na.rm = TRUE) / ( sqrt( sum( x1^2 , na.rm = TRUE)) * sqrt( sum( x2^2 , na.rm = TRUE) ) )
  # TODO: to address div by zero if one vector is all 0, could use a small value 
  # "epsilon" (1e-8) to denomintor if zero

  # if one of the sites doesn't have any individuals/species in it
  # i.e. the row sum is 0
  
  # TODO: ensure all values between 0 and ?1?!
  if (cosine < 0){ cosine <- 0}
  if (cosine > 1){ cosine <- 1}
  
  
  # # TODO: think I should scrap the below adjustments...
  # if ((sum(x1, na.rm = TRUE) == 0) |
  #     (sum(x2, na.rm = TRUE) == 0)) {
  #   # then the similarity between sites should be 0
  #   cosine <- 0
  #   # if both sites have no individuals/species
  # } else if (((sum(x1, na.rm = TRUE) == 0) &
  #             (sum(x2, na.rm = TRUE) == 0)) ){
  #   # then class the similarity as NA
  #   cosine <- NA
  # }
  cosine
}


#' Calculate horn overlap metric - by default, a similarity
#' 
#' @description 
#' This function calculates Horn's overlap (1966) between abundances of 
#' species at two sites
#' @details 
#' Expects two vectors, one each for the abundance/occurence of species at two sites
#' assuming that values at the same index (e.g., x[i], y[i]) correspond to the 
#' same species/taxonomic entity
#' @param x A [`numeric`] vector for abundance/occurence of species at site 1
#' @param y A [`numeric`] vector for abundance/occurence of species at site 2
#' @keywords database, processing
#' @examples 
#' \dontrun{
#'   horn_overlap <- calc_horn_overlap(c(0,1,1), c(1,0,1))
#' } 
#' @noRd
calc_horn_overlap <- function(x, y){
  # xp <- x / sum(x, na.rm = T)
  # yp <- y / sum(y, na.rm = T)
  num1 <- sum( (x + y) * log(x + y), na.rm = T )
  num2 <- sum(x * log(x), na.rm = T)
  num3 <- sum(y * log(y), na.rm = T)
  
  # print(paste("Num1:", num1))
  # print(paste("Num2:", num2))
  # print(paste("Num3:", num3))
  
  num <- num1 - num2 - num3
  
  x_sum <- sum(x, na.rm = T)
  y_sum <- sum(y, na.rm = T)
  xy_sum <- x_sum + y_sum
  
  denom1 <- xy_sum * log(xy_sum)
  denom2 <- x_sum * log(x_sum)
  denom3 <- y_sum * log(y_sum)
  
  # print(paste("Denom1:", denom1))
  # print(paste("Denom2:", denom2))
  # print(paste("Denom3:", denom3))
  
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


#' Calculate non-vegan-based dissimilarity metrics
#' 
#' @description 
#' This function calculates non-vegan-based dissimilarity metrics between sites
#' @details 
#' Expects a dataframe/tibble with species abundance records at two sites 
#' - sites as rows, species as columns
#' @param spp_site_data A dataframe/tibble object.
#' @param metric A [`character`] - one of "hypervolume", "cosine", "theta",
#' "baselga_bray_balanced", "baselga_bray_gradient", "baselga_ruzicka_balanced", 
#' "baselga_ruzicka_gradient".
#' @param binary A [`logical`] - should distance be based on presence/absence
#' @param rel_abund A [`logical`] - whether relative abundance values should 
#' be used; applicable for hypervolume and baselga_ metrics. Cosine and theta are 
#' unaffected by relative abundance v raw abundance.
#' @param rescale A [`logical`] - should distance be rescaled relative to
#' theoretical max (if applicable and not 1).
#' @param div_metric_type A [`character`] describing the type of biodiversity records
#' at the site of interest - one of "Occurrence" or "Abundance". Set by hib_cacl_div_between.
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
                            rescale,
                            div_metric_type
                            ){
  
  assertthat::assert_that(
    is.character(metric),
    is.logical(binary) || is.null(binary),
    is.logical(rel_abund) || is.null(rel_abund),
    is.logical(rescale) || is.null(rescale),
    is.character(div_metric_type)
    # assertthat::has_name(spp_site_data, "SSBS"),
    # assertthat::has_name(spp_site_data, "Taxon_name_entered"),
    # assertthat::has_name(spp_site_data, "Measurement")
  )
  # print(use_relative_abundance)
  metric <- match.arg(metric, 
                      c("hypervolume", "cosine", "theta",
                        "baselga_bray_balanced", "baselga_bray_gradient", 
                        "baselga_ruzicka_balanced", "baselga_ruzicka_gradient",
                        "horn_dissim"),
                      # "baselga_bray","baselga_ruzicka",
                      # "baselga_jaccard", "baselga_sorensen") 
                      several.ok = FALSE)
  # TODO: challenge of single spp communities?? - this would affect metric calcs..??

  # TODO:
  # if!binary, ensure only abundance data is used...
  if ((binary == F) & (div_metric_type == "Occurrence")){
    return(NULL)
  }
  # Convert to relative abundance if needed
  if (rel_abund){
    if (metric %in% c("hypervolume", 
                        "baselga_bray_balanced", "baselga_bray_gradient", 
                        "baselga_ruzicka_balanced", "baselga_ruzicka_gradient",
                        "horn_dissim")){
      row_abund <- rowSums(spp_site_data, na.rm = TRUE)
      row_abund[row_abund == 0] <- 1
      spp_site_data <- spp_site_data/row_abund
    }
  }
  # Use species occurence if specified
  # if!binary, enruse only abundance data is used...
  if (binary){
    spp_site_data[spp_site_data > 0] <- 1
  }

  # Add in hypervolume approach
  if (metric == "hypervolume"){
    
    # transpose site_spp_data for use in convhulln (M rows of N dimensions)
    spp_site_data <- t(spp_site_data)
    
    # add origin to spp_site_data
    spp_site_data <- rbind(spp_site_data,
                           rep(0, ncol(spp_site_data)))

    # Should an error default to 0
    # 0 should be vol if either site has 0 total abund
    # And technically if both comms are 0..
    vol <- tryCatch(
      {
        geometry::convhulln(spp_site_data, 
                            output.options = TRUE)$vol
      },
      error = function(e) {
        # NA # ??
        0
        # Wouldn't NA be a better error value
      }
    )
    
    if (rescale & rel_abund){
      vol <- vol / 0.5
    }
    
    out <- data.frame("Hypervolume" = vol)
  }
  
  
  else if (metric %in% c("cosine", "theta")){
    
    x <- as.matrix(spp_site_data)
    # N rows == 2
    x1 <- x[1,]
    x2 <- x[2,]

    cosine_sim <- calc_cosine_sim(x1, x2)

    if (metric == "cosine"){
      # TODO: ensure all rescaled values between 0 and 1!
      out <- data.frame("Cosine" = 1 - cosine_sim)
    }
    if (metric == "theta"){
      theta <- acos(cosine_sim)
      if (theta < 0){ 
        theta <- 0
      }
      if (rescale){
        theta <- theta / (pi / 2)
        # Ensure all rescaled values between 0 and 1!
        if (theta > 1){ theta <- 1}
      }
      # Else ensure 0 - pi/2
      if (theta > (pi/2)){ 
        theta <- (pi/2)
      }
      out <- data.frame("Theta" = theta)
    }
  }

  else if (metric == "horn_dissim"){
    horn_overlap <- calc_horn_overlap(x1, x2)
    out <- data.frame("Horn_dissim" = 1 - horn_overlap)
  }
  
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
#' Expects a [`HibDatabase`] object.
#' @param obj A [`HibDatabase`] object.
#' @param metrics A [`character`] vector of dissimilarity metrics. 
#' Metrics available via \code{vegan::vegdist}: "manhattan", "euclidean", "canberra", "clark", 
#' "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", 
#' "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", 
#' "hellinger", "aitchison", "robust.aitchison".
#' Other available metrics: "hypervolume", "cosine", "theta", "baselga_bray_balanced", 
#' "baselga_bray_gradient", "baselga_ruzicka_balanced", "baselga_ruzicka_gradient",
#' "horn_dissim".
#' Metrics can be provided in upper or lower case - internal processing is done 
#' using lower-case - and as as unambiguous abbreviations.
#' Note: "horn" refers to the "Morisita-Horn" index in the vegan package. 
#' "horn_dissim" refers to (1 - Horn overlap), calculated endogenously within this package.
#' @param binary A [`logical`] vector, should metric be calculated using occurence
#' values (if TRUE) or abundance (if FALSE, default).
#' @param rel_abund A [`logical`] vector determining whether relative abundance 
#' should be used to calculate the dissimilarity metric, only applicable for some metrics, 
#' e.g., euclidean, hypervolume - Note: many metrics include rescaling to 'relative abundances'
#' as part of the calculation. 
#' Default values is TRUE.
#' @param rescale A [`logical`] vector determining if dissimilarity measure should 
#' be rescaled relative to its maximum value - applies to metrics that are not 
#' on a 0-1 scale by default - hypervolume, euclidean (if relative abundance), 
#' chord, hellinger
#' Default values is "Primary vegetation_Minimal use" - the most pristine LUI within
#' the PREDICTS database.
#' @param baseline_lui A [`character`] or NULL, indicating the baseline LUI against
#' which compositional dissimilarity is calculated per study. If NULL, compositional 
#' dissimilarity is calculated between all sites within each study.
#' Default values is 'Primary vegetation_Minimal use'.
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
                                 binary = FALSE,
                                 rel_abund = TRUE,
                                 rescale = TRUE,
                                 baseline_lui = "Primary vegetation_Minimal use"
                                 ){

  # TODO:
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
    
    # assertthat::has_name(obj$species, "Predominant_land_use"),
    # assertthat::has_name(obj$species, "Use_intensity"),
    # assertthat::has_name(obj$species, "LUI"),
    # TODO: if not present, they have to be in obj$sites, hib_create...
    # Can left join in - LUI is essential
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
                       "horn_dissim")
  # "baselga_bray","baselga_ruzicka",
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
  # print(metrics_spec_df)
  # print(metrics_)
  matchable_metrics_idx <- !is.na(pmatch(metrics_spec_df$metric,
                                         metrics_))
  # print(matchable_metrics_idx)
  # subset to matchable metrics
  metrics_spec_df <- metrics_spec_df[matchable_metrics_idx,]
  # print(metrics_spec_df)
  # run match.arg on matchable metrics
  metrics_spec_df$metric <- match.arg(metrics_spec_df$metric, 
                                       metrics_,
                                       several.ok = T)
  
  metrics_spec_df <- unique(metrics_spec_df)
  # print(metrics_spec_df)
  # message(paste(metrics_spec_df$metric, collapse = ","))

  # %%%
  # if metric is already in hib$diversity_between, ignore
  # TODO: check recalculate...
  if (!is.null(obj$diversity_between)){

    # Need to regen altGower
    metrics_tc <- metrics_spec_df$metric
    # print(metrics_tc)
    # metrics_tc[metrics_tc == "altgower"] <- "altGower"
    metrics_tc <- tools::toTitleCase(metrics_tc)
    # print(metrics_tc)
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
    # if nrow(metrics_spec_df) == 0
    # print(nrow(metrics_spec_df))
    if (nrow(metrics_spec_df) == 0){
        # message(crayon::red(paste0("\nError: no uncalculated biodiversity metrics.\n\n",
        #                        "Check that the specified biodiversity responses ('", baseline_lui,
        #                        "') is present ...\n"))
        #   )
        stop("No new biodiversity metrics to calculate.", call. = F)
    }
  }
  # %%%
  
  message(crayon::green("\nPreparing site-contrasts to calculate compositional dissimilarity ...\n"))
  

  
  # Get site comparisons

  # TODO: if all metrics are binary 
  # can also use 'Occurrence' data...
  # Diversity_metric_type == 'Occurrence'
  # if (all(metrics_spec_df$binary)){}

  # if (!all(metrics_spec_df$binary)){}
  # Or could include occurence, and filter it out in subsequent functions if !binary...

  spp_tmp <- obj$species |>

  # TODO:
    # dplyr::filter(#!is.na(Predominant_land_use) & # ?? - drop this
    #               Diversity_metric_type == "Abundance") |>

    dplyr::filter(#!is.na(Predominant_land_use) & # ?? - drop this
                  Diversity_metric_type %in% c("Abundance",
                                               "Occurrence")) |>
    
    dplyr::group_by(SS) |>
    dplyr::mutate(n_sample_effort = dplyr::n_distinct(Sampling_effort),
                  # Taxon_name_entered v best_guess_binomial ??
                  n_species = dplyr::n_distinct(Taxon_name_entered)) |>
    dplyr::ungroup() |>
    dplyr::filter(n_sample_effort == 1 &
                    n_species > 1)
  
  # print(metrics_spec_df)
  if (all(metrics_spec_df$binary == FALSE)){
    spp_tmp <- spp_tmp |>
      dplyr::filter(Diversity_metric_type == "Abundance")

    # print(nrow(spp_tmp))
    if (nrow(spp_tmp) == 0){
      stop("\nError: All metric specifications require Abundance records (binary = FALSE), but none are available.", call. = F)
    }

  }

  
  # print(unique(spp_tmp$Diversity_metric_type))
  if (!is.null(baseline_lui)){
    
    # Better appr is to first get site data..., (subset to spp_tmp_SSBS) 
    # then check for baseline LUI...

    # %%%%
    site_tmp <- obj$sites |>
      dplyr::filter(SSBS %in% unique(spp_tmp$SSBS),
                    !is.na(LUI))
    # TODO: filter for !is.na(LUI) here!

    if (nrow(site_tmp) == 0){
      stop("\nError: no site comparisons possible. Check your data and the metric specifications.", call. = F)
    }

    assertthat::assert_that(baseline_lui %in% unique(site_tmp$LUI),
                            msg = "baseline_lui is not present in the site data that can be used for comparisons. \n      Consider changing your data and/or the metric specifications.")

    site_tmp$base_lui <- site_tmp$LUI == baseline_lui
    
    site_tmp <- site_tmp |>
      dplyr::group_by(SS) |>
      dplyr::mutate(n_base_lui = sum(base_lui)) |>
      dplyr::filter(n_base_lui > 0) |>
      dplyr::ungroup() |>
      droplevels()
    # %%%%

    # This is an issue if LUI updates have only been added to the site, not spp data!!
    
    # spp_tmp$base_lui <- spp_tmp$LUI == baseline_lui
    
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
  
  # print(unique(site_tmp$Diversity_metric_type))

  n_spp <- spp_tmp |>
    dplyr::select(Best_guess_binomial) |>
    dplyr::n_distinct()
  n_spp_abund <- spp_tmp |>
    dplyr::filter(Diversity_metric_type == "Abundance") |>
    dplyr::select(Best_guess_binomial) |>
    dplyr::n_distinct()

  # TODO: n_spp_abund, n_spp_abund_occ

  # TODO: n taxa ?
  # TODO: drop this info (spp/taxa) here ? - not really useful??  

  # TODO: set any NA measures to 0 for safety
  if (sum(is.na(spp_tmp$Measurement)) > 0){
    message(crayon::yellow(paste0("\nWarning: NAs detected in 'Measurement'.  Setting NAs to zeros.")))
    spp_tmp$Measurement[is.na(spp_tmp$Measurement)] <- 0
  }

  assertthat::assert_that(
    sum(is.na(spp_tmp$Measurement)) == 0,
    is.numeric(spp_tmp$Measurement)
  )
  # spp_tmp$Measurement[is.na(spp_tmp$Measurement)] <- 0

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
  
  # if nrow(site_comp_df == 0) stop...
  if (nrow(site_comp_df) == 0){
    message(crayon::red(paste0("\nError: no site comparisons possible...\n\n",
                               "Maybe check that the specified baseline LUI ('", baseline_lui,
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

  # print(head(site_comp_df))

  site_comp_df <- site_comp_df |>
    dplyr::select(-Diversity_metric_type_s2) |>
    dplyr::rename(Diversity_metric_type = Diversity_metric_type_s1)



  # Some checks of metric specs and diversity data types available
  # print(unique(site_comp_df$Diversity_metric_type))
  un_div <- unique(site_comp_df$Diversity_metric_type) |> as.character()
  # TODO: if all binary are false, / 
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
                                   
                                   # print( spp_data |>
                                   #   # filter out the SS that matches the SS we're interested in
                                   #   dplyr::filter(SS == ss) |>
                                   #   # pull out the site name, species name and abundance information
                                   #   dplyr::select(SSBS, Taxon_name_entered, Measurement) |>
                                   #   # pivot so that each column is a species and each row is a site
                                   #   tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement,
                                   #                      values_fill = 0) )
                                     # TODO: fix values_fill bug...
                                     # tidyr::pivot_wider(names_from = Taxon_name_entered, 
                                     #                    values_from = Measurement) )


                                   spp_site_data <- spp_data |>
                                     # filter out the SS that matches the SS we're interested in
                                     dplyr::filter(SS == ss) |>
                                     # pull out the site name, species name and abundance information
                                     # TODO: could consider Diversity_metric_type here...!
                                     # pass to hib_calc_veg...
                                     dplyr::select(SSBS, Taxon_name_entered, Measurement) |>                                  
                                     # pivot so that each column is a species and each row is a site
                                     tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement,
                                                        values_fill = 0) |>
                                     
                                     # set the rownames to the SSBS and then remove that column
                                     tibble::column_to_rownames("SSBS")
                                   
                                   div_metr_type <- comp_df |>
                                     dplyr::filter(SS == ss) |>
                                     dplyr::pull(Diversity_metric_type) |>
                                     unique() |> as.character()

                                   # Get Diversity metric_type of each SSBS
                                   # div_tmp <- spp_data |>
                                   #   dplyr::filter(SS == ss) |>
                                   #   dplyr::select(SSBS, Diversity_metric_type) |>
                                   #   unique() 
                                   
                                   # if ("Occurrence" %in% div_tmp$Diversity_metric_type){
                                   #  print(div_tmp)
                                   #  print(all(div_tmp$Diversity_metric_type == "Occurrence"))
                                   # }

                                   # div_metr_type <- unique(div_tmp$Diversity_metric_type)
                                   
                                   assertthat::assert_that(
                                     length(div_metr_type) == 1
                                   )

                                   # print(div_metr_type)

                                   # TODO: check at least 2 species are present (i.e 2 cols or more)
                                   # TODO: return max_dist value - if NA, set rescale to F...
                                   # Depending on metric_vec values, calc metrics of interest
                                   # %%%
                                   out_df <- purrr::pmap(.l = metric_df,
                                                            .f = function(metric,
                                                                          binary,
                                                                          rel_abund,
                                                                          rescale,
                                                                          spp_data = spp_site_data,
                                                                          div_metric_type = div_metr_type){
                                                              # print(metric)
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
  
  # do other metrics...
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
                                  # TODO: could consider Diversity_metric_type here...!
                                     # pass to hib_calc_nonveg...
                                  tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement,
                                                     values_fill = 0) |>
                                  # set the rownames to the SSBS and then remove that column
                                  tibble::column_to_rownames("SSBS")
                                
                                # Get Diversity metric_type of each SSBS
                                # spp_data |>
                                #   dplyr::filter(SSBS %in% c(s1, s2)) |>
                                #   dplyr::select(SSBS, Diversity_metric_type) |>
                                #   unique() |> print()


                                ss <- spp_data |>
                                  # filter out the SSBS that matches the pair of sites we're interested in
                                  dplyr::filter(SSBS %in% c(s1, s2)) |>
                                  dplyr::pull(SS) |>
                                  unique()
                                
                                out_df <- data.frame(SSBS1 = s1, 
                                                     SSBS2 = s2,
                                                     SS    = ss)
                                
                                # TODO: non-binary metrics..!!!
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
  
  # convert to long tibble...
  # add processing flags per metric - 
  
  # TODO: add binary to veg df...
  # If binary == False, set any metric to NA where Diversity metric type == Occurrence



  # Calc is_empty flags once and add to each sub df?
  empty_sites1 <- spp_tmp |>
    dplyr::group_by(SSBS) |>
    dplyr::select(SSBS, Taxon_name_entered, Measurement) |>
    # pivot so that each column is a species and each row is a site
    tidyr::pivot_wider(names_from = Taxon_name_entered, values_from = Measurement,
                       values_fill = 0) |>
    tibble::column_to_rownames("SSBS") #|>
    # dplyr::mutate(SSBS1_is_empty = (rowSums(.) == 0)) |>
    # tibble::rownames_to_column("SSBS1") |>
    # dplyr::select(SSBS1, SSBS1_is_empty)
    # %%%%%

  # TODO: ignore NAs 
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
  # print(head(veg_df))
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
    dplyr::mutate(baseline_lui = baseline_lui)
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
  
  if (rel_abund == T){
    veg_max_df$Theoretical_max[veg_max_df$Metric == "Euclidean"] <- sqrt(2)
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
  # hypervol = NA, or 0.5 if rel_abund
  # cosine = 0-1 (set 1 as max, convert to dissim)
  # theta  = pi/2
  # all others, 0-1 - check with baselga...
  non_veg_max_df <- data.frame(Metric = tools::toTitleCase(c("hypervolume", "cosine", "theta",
                                                             "baselga_bray_balanced", "baselga_bray_gradient", 
                                                             "baselga_ruzicka_balanced", "baselga_ruzicka_gradient",
                                                             "horn_dissim")),
                              Theoretical_max = c(NA, 1, pi/2,
                                                 1, 1,
                                                 1, 1,
                                                 1))
  if (rel_abund == T){
    non_veg_max_df$Theoretical_max[non_veg_max_df$Metric == "Hypervolume"] <- 0.5
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

