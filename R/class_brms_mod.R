#' Template objects for brms model 
#' 
#' @description 
#' This is the template object for a brms model class.
#' It contains all the functions inherent to the model class. Any functions in 
#' the package that operate on or interact with the model(summarising, plottting) 
#' need to be coded in here or provided separately.
#' @keywords internal
#' @import R6
#' @noRd
hib_brms_model <- R6Class("hib_brms_model",
                          # Public data objects and functions
                          public = list(
                            version = NULL,
                            last_update = NULL,
                            # paths = NULL, # ??
                            data = NULL, # - restrict data to that used in model fitting... + version?
                            model = NULL,
                            extra_args = NULL,
                            fixedef_coefficients = NULL,
                            metrics = NULL,
                            processing_flags = NULL,
                            # TODO: what else to include
                            # 
                            # lme4 output is pretty comprehensive, this is just a wrapper...
                            # need optimal ranef struct, if found
                            # need model ranks, if dredged..
                            
                            # Executed when created
                            initialize = function(db_obj, # ??
                                                  mod_df,
                                                  brms_obj,
                                                  fixef_df,
                                                  metr_ls,
                                                  n_spp,
                                                  extra_args
                                                  # set_up_specs
                                                  # spp df?
                                                  # ???? - should this do the model fitting
                                                  # or should that be in hib_fit_lme4
                                                  # and this just return the resultant obj
                                                  # facilitating summary, printing...
                            ){
                              # - some assertion checks
                              assertthat::assert_that(
                                # Check database obj
                                # db_obj$checkself(),
                                inherits(db_obj, "Database") || is.list(db_obj),
                                # Check model data
                                is.data.frame(mod_df),
                                # Check model
                                class(brms_obj) == "brmsfit",
                                # Check coefficient df and metrics
                                is.data.frame(fixef_df),
                                is.list(metr_ls),
                                is.numeric(n_spp),
                                is.list(extra_args)
                              )
                              
                              # ------ #
                              message(crayon::green("Initialising fitted hib model object ...") )
                              
                              # %%%%%%%%%%%%%%%
                              
                              # Save inputs in self
                              # Set mod_df as data
                              self$data <- mod_df
                              # Set model as brms_obj
                              self$model <- brms_obj

                              self$extra_args <- extra_args
                              
                              # Set fixed_coefficients as fixef_df
                              self$fixedef_coefficients <- fixef_df
                              # Set metrics as metr_ls
                              self$metrics <- metr_ls
                              
                              # Incl processing flags from hib_db
                              self$processing_flags <- db_obj$processing_flags
                              # TODO: what about additional data (data2), e.g. covariance structures
                              
                              # Calc n sites/species used in  model - need to calc spp with div metrics!!
                              # N_sites <- dplyr::n_distinct(mod_df$SSBS)
                              
                              if ("SSBS" %in% colnames(mod_df)){
                                self$processing_flags$N_sites <- dplyr::n_distinct(mod_df$SSBS)
                              } else if ("SSBS1" %in% colnames(mod_df)){
                                # print("SSBS1 found")
                                self$processing_flags$N_sites <- length(unique(c(mod_df$SSBS1, mod_df$SSBS2)))
                              }
                              
                              self$processing_flags$N_spp <- n_spp
                              
                              if (inherits(db_obj, "Database")){
                                self$version <- db_obj$version
                              }
                              if (is.list(db_obj)){
                                self$version <- lapply(db_obj, function(x){x$version})
                              }
                              
                              self$last_update <- lubridate::today()
                              
                              # %%%%%%%%%%%%%%%
                           
                            },
                            # Functions
                            
                            # Base printout function
                            print = function(){
                              # Padding numbers function, from:
                              # https://github.com/jthomasmock/gtExtras/blob/HEAD/R/pad_fn.R
                              pad_fn <- function(x, nsmall = 2, pad0) {
                                # round and format values as text with specific number of decimals
                                round_x <- round(x, digits = nsmall)
                                fmt_x <- format(round_x, nsmall = nsmall)
                                
                                # calc number of trailing zeros
                                zero_len <- nchar(fmt_x) - nchar(sub("0*$", "", fmt_x))
                                
                                # create string of zero
                                if (pad0) {
                                  rep_zero <- strrep(0, zero_len)
                                } else {
                                  rep_zero <- ""
                                }
                                
                                # remove trailing zeros by position
                                fmt_out <- substr(fmt_x, 1, nchar(fmt_x) - zero_len)
                                
                                # add the non-breaking spaces to the formatted values
                                filled_out <- paste0(fmt_out, rep_zero)
                                
                                filled_out
                              }
                              # N spp/sites...
                              # print(summary(self$model))
                              stats <- list(number_sites = self$processing_flags$N_sites,
                                            # number_lu = length( unique(self$data$Predominant_land_use) ),
                                            # number_lui = length(unique(self$data$Use_intensity)),
                                            number_species = paste(self$processing_flags$N_spp, collapse = ", ") ) #,
                              # filters_applied = paste(self$processing_flags$filtering_steps, collapse = "; "),
                              # div_metrics = paste(sort(unique(self$diversity$Metric)), collapse = ", ") )  
                              
                              message(
                                paste0("Hibiiscus fitted brms model\n",
                                       "(Last updated: ", crayon::green( self$last_update ), ")",
                                       "\n",
                                       "\t Number of sites:               ", crayon::green(stats$number_sites), "\n",
                                       # "\t Land use (class x intensity):  ", crayon::green(stats$number_lu), " x ",
                                       # crayon::green(stats$number_lui), "\n",
                                       "\t Number of species:             ", crayon::green(stats$number_species), "\n",
                                       "\n\n",
                                       "--------------------\n\n",
                                       "Summary of fitted model:\n"
                                       # "Fixed effects:"
                                )
                              )
                              print(summary(self$model))
                              
                              # Get loo/waic...
                              metr <- c(self$metrics$loo$estimates["looic","Estimate"],
                                        self$metrics$waic$estimates["waic","Estimate"])
                              metr <- pad_fn(metr, 3, F)
                              names(metr) <- c("looic", "waic")
                              
                              message(paste0("\n\n",
                                             # "--------------------\n",
                                             # "Random effects:\n",
                                             # ranef_chr, "\n\n",
                                             "Performance metrics:\n",
                                             "\tLOO IC \t ", metr["looic"], "\n",
                                             "\tWAIC   \t ", metr["waic"], "\n"))
                              
                            },
                            # Save object
                            save_hib = function(output_path = NULL){
                              assertthat::assert_that(
                                is.null(output_path) || is.character(output_path)
                              )
                              
                              if(!assertthat::has_extension(output_path, "rds")){
                                output_path <- paste0(output_path, ".rds")
                              }
                              
                              if (file.exists(output_path)){
                                message(crayon::yellow(paste0("Warning: '", output_path, "' already exists.",
                                                              "\n",
                                                              "This file will be updated.\n")))
                              }
                              
                              message(crayon::green(paste0("Saving hib model to: '", output_path, "' ...\n")) )
                              
                              saveRDS(self, output_path)
                              assertthat::assert_that(file.exists(output_path),
                                                      msg = "Writing a copy of the model failed...")
                              invisible()
                            },
                            
                            plot = function(type = "fixed-effects"){
                              # sjPlot::plot_model(self$model, 
                              #                    show.intercept = T)
                              
                              assertthat::assert_that(is.character(type))
                
                              # Plot options
                              type = match.arg(type,
                                               c("fixed-effects",
                                                 "random-effects"))
                              
                              basic_thm <- ggplot2::theme(axis.text = ggplot2::element_text(size = 16),
                                                 axis.title = ggplot2::element_text(size = 20),
                                                 plot.title = ggplot2::element_text(size = 20),
                                                 plot.subtitle = ggplot2::element_text(size = 18),
                                                 strip.text.x = ggplot2::element_text(size = 18),
                                                 strip.text.y = ggplot2::element_text(size = 18),
                                                 strip.background = ggplot2::element_blank(),
                                                 legend.text = ggplot2::element_text(size = 16),
                                                 legend.title = ggplot2::element_text(size = 18),
                                                 legend.text.align = 0) 
                              
                              # 1. fixed effects - summary
                              if (type == "fixed-effects"){
                                fixef_df <- self$fixedef_coefficients
                                
                                # Set a baseline- Intercept
                                fixef_df$Coefficient_name <- factor(fixef_df$Coefficient_name)
                                if ("Intercept" %in% fixef_df$Coefficient_name){
                                  fixef_df$Coefficient <- relevel(fixef_df$Coefficient_name,
                                                                  ref = "Intercept")
                                }
                                
                                # %%%
                                plt <- ggplot2::ggplot(fixef_df) +
                                  ggplot2::geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                                  ggplot2::geom_errorbar(ggplot2::aes(x = Coefficient_name, 
                                                    ymin = loCI,
                                                    ymax = hiCI),
                                                width = 0) +
                                  ggplot2::geom_point(ggplot2::aes(x = Coefficient_name,
                                                 y = Coefficient_value)) +
                                  
                                  ggplot2::labs(y = "Estimate",
                                       x = ggplot2::element_blank(),
                                       title = "Fixed effect estimates") +
                                  ggplot2::theme_bw() +
                                  basic_thm +
                                  ggplot2::coord_flip()
                                
                              }
                              #    fixed effects - posteriors
                              # 2. random iid
                              if (type == "random-effects"){
                                # TODO: ranef(m) produces a list of 3D arrays
                                # dim 1 = factor levels
                                # dim 2 = summary stats
                                # dim 3 = group-level effects
                                
                                # if random slopes, then x is 3d, not 2d
                                # coef nm, (intercept or LU) also in colnames
                                rand_df <- purrr::map(brms::ranef(self$model), 
                                           function(x){
                                             
                                             # melt - then colnames are coef_level, coef_stat, coef_term
                                             tmp <- reshape2::melt(x) %>%
                                               # widen across stat
                                               tidyr::pivot_wider(names_from = "Var2",
                                                                  values_from = "value") %>%
                                               dplyr::rename(Coef_level = Var1,
                                                             Coef_term = Var3,
                                                             Coefficient_value = Estimate,
                                                             Coefficient_error = Est.Error,
                                                             loCI = Q2.5,
                                                             hiCI = Q97.5)
                                             # return
                                             # tmp <- as.data.frame(x)
                                             # colnames(tmp) <- c()
                                             tmp
                                             }) %>%
                                  dplyr::bind_rows(.id = "Coefficient_group") %>%
                                  # tibble::rownames_to_column("Coefficient_name") %>%
                                  dplyr::group_by(Coefficient_group, Coef_term) %>%
                                  dplyr::arrange(Coefficient_value, .by_group = T) %>%
                                  dplyr::mutate(ID_ord = 1:dplyr::n()) %>%
                                  dplyr::ungroup() %>%
                                  dplyr::mutate(sig_col = dplyr::case_when(loCI > 0 ~ "+ve",
                                                                    hiCI < 0 ~ "-ve",
                                                                    T ~ "n.s."))
                                
                                
                                # colour as +ve or -ve 
                                # plot - no names...
                                plt <- ggplot2::ggplot(rand_df) +
                                  ggplot2::geom_hline(yintercept = 0, lty = 2, colour = "darkgrey") +
                                  ggplot2::geom_errorbar(ggplot2::aes(x = ID_ord, 
                                                    ymin = loCI,
                                                    ymax = hiCI,
                                                    colour = sig_col),
                                                show.legend = F,
                                                alpha = 0.5,
                                                width = 0) +
                                  ggplot2::geom_point(ggplot2::aes(x = ID_ord,
                                                 y = Coefficient_value,
                                                 colour = sig_col),
                                             show.legend = F,
                                             alpha = 0.7) +
                                  ggplot2::labs(y = "Estimate",
                                       x = ggplot2::element_blank(),
                                       title = "Random effect estimates") +
                                  ggplot2::scale_color_manual(values = c("+ve" = "dodgerblue3",
                                                                "-ve" = "firebrick3",
                                                                "n.s."= "grey30")) +
                                  ggplot2::facet_grid(Coef_term~Coefficient_group, scales = "free_x") +
                                  ggplot2::theme_bw() +
                                  basic_thm +
                                  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                        aspect.ratio = 1)

                              }
                              (plt)
                            }
                          ),
                          private = list()
)