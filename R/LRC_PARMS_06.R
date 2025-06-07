#' @title Fit a Light Response Curve (LRC) by an index to get a parameter file
#'
#' @description
#' This function uses the equation:
#' \deqn{\text{NEE} \sim a_x \left( 1 - \exp \left( a_1 \cdot \left( \text{PAR} - I_{\text{comp}} \right) \right) \right)}
#'
#' Where NEE is daytime net ecosystem exchange (\eqn{\mu}mol CO2 m-2 s-1),
#' \eqn{a_x} is the maximum photosynthesis rate (\eqn{\mu}mol CO2 m-2 s-1),
#' \eqn{a_1} is a coefficient,
#' and \eqn{I_{\text{comp}}} is the compensation light level with the same unit of PAR.
#'
#' The equation requires photosynthetically active radiation, PAR, in \eqn{\mu}mol m-2 s-1 and net ecosystem exchange, NEE, in \eqn{\mu}mol m-2 s-1.
#'
#' @details
#' Model parameters are fit using the R package `brms`.
#'
#' Rhat (Potential Scale Reduction Factor):
#' Indicates how well the different Markov chains in your analysis have converged to the same posterior distribution.
#' Ideally, Rhat should be close to 1 for all parameters.
#' A high Rhat value suggests potential convergence issues and the need to run the chains longer.
#'
#' Bulk ESS (Effective Sample Size - Bulk):
#' Estimates the effective number of independent samples from the central part of the posterior distribution.
#'
#' Tail ESS (Effective Sample Size - Tail):
#' Estimates the effective number of independent samples from the tails of the posterior distribution.
#' Important for assessing the reliability of quantile estimates (e.g., 95% confidence intervals).
#'
#' Key points to remember:
#' Aim for Rhat close to 1 and high values for both Bulk ESS and Tail ESS.
#'
#' @param data.frame (dataframe) A dataframe that contains net ecosystem exchange (NEE), an index, and photosynthetically active radiation (PAR).
#' @param iterations (numeric) The number of iterations to run `brms::brm()`.
#' @param priors.lrc (brmsprior dataframe) The priors for `brms::brm()` to use.
#' Default priors are as follows:
#' ```
#' brms::prior("normal(-5, 5)", nlpar = "ax", lb = -30, ub = 0) +
#' brms::prior("normal(0.1, 1)", nlpar = "alpha", lb = 0, ub = 1) +
#' brms::prior("normal(300, 100)", nlpar = "Icomp", lb = 0, ub = 500)
#' ```
#' @param idx.colname (character) The name of the column containing the index.
#' @param NEE.colname (character) The name of the column containing NEE.
#' @param PAR.colname (character) The name of the column containing PAR.
#'
#' @returns (dataframe) Dataframe of parameter values by the index used to fit them.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Set the working directory to the location of the sampledata file: AMF_US-Skr_BASE_HH_2-5_Formatted.csv
#'setwd('sampledata')
#'# Import flux tower data
#'tower.data <- read.csv('AMF_US-Skr_BASE_HH_2-5_Formatted.csv')
#'# Fit curve parameters for each YearMon:
#'Example_LRC_PARMS_06 <-LRC_PARMS_06(data.frame = tower.data,
#'                                    iterations = 1000,
#'                                    priors.lrc = brms::prior("normal(-5, 5)", nlpar = "ax", lb = -30, ub = 0) +
#'                                      brms::prior("normal(0.1, 1)", nlpar = "alpha", lb = 0, ub = 1) +
#'                                      brms::prior("normal(300, 100)", nlpar = "Icomp", lb = 0, ub = 500),
#'                                    idx.colname = 'YearMon',
#'                                    NEE.colname = 'NEE_PI',
#'                                    PAR.colname = 'SW_IN')
#'
LRC_PARMS_06 <- function(data.frame = NULL,
                         iterations = NULL,
                         priors.lrc = brms::prior("normal(-5, 5)", nlpar = "ax", lb = -30, ub = 0) +
                           brms::prior("normal(0.1, 1)", nlpar = "alpha", lb = 0, ub = 1) +
                           brms::prior("normal(300, 100)", nlpar = "Icomp", lb = 0, ub = 500),
                         idx.colname = NULL,
                         NEE.colname = NULL,
                         PAR.colname = NULL){

  # Squelch visible bindings note
  nee <- idx <- PAR <- NULL

  data.frame$nee <- data.frame[,NEE.colname]
  data.frame$idx <- data.frame[,idx.colname]
  data.frame$PAR <- data.frame[,PAR.colname]

  df <- data.frame %>% dplyr::select(idx,
                                     nee,
                                     PAR)

  equation <- nee ~ ax * (1 - exp(alpha * (PAR - Icomp)))


  if(c("idx") %in% base::names(df)) {
    base::print("GREAT JOB! your dataframe contains idx")
  } else {base::print("The dataframe must include: idx, nee, and PAR")}

  if(c("nee") %in% base::names(df)) {
    base::print("YIPEE! your dataframe contains nee")
  } else {base::print("The dataframe must include: idx, nee, and PAR")}

  if(c("PAR") %in% base::names(df)) {
    base::print("Hooray! your dataframe contains PAR")
  } else {base::print("The dataframe must include: idx, nee, and PAR")}

  # PARM Dataframe:
  parms <- base::data.frame(idx = base::as.character(),
                            ax.mean = base::as.numeric(),
                            ax.se = base::as.numeric(),
                            ax.Bulk_ESS = base::as.numeric(),
                            ax.Tail_ESS = base::as.numeric(),
                            ax.Rhat = base::as.numeric(),

                            alpha.mean = base::as.numeric(),
                            alpha.se = base::as.numeric(),
                            alpha.Bulk_ESS = base::as.numeric(),
                            alpha.Tail_ESS = base::as.numeric(),
                            alpha.Rhat = base::as.numeric(),

                            Icomp.mean = base::as.numeric(),
                            Icomp.se = base::as.numeric(),
                            Icomp.Bulk_ESS = base::as.numeric(),
                            Icomp.Tail_ESS = base::as.numeric(),
                            Icomp.Rhat = base::as.numeric(),
                            samples = base::as.numeric())

  base::message(" Your dataframe looks good and you are now ready to start fitting models")

  for (i in base::unique(df$idx)){
    base::print(i)

    # Subset the file:
    df.sub <- df %>% dplyr::filter(idx == i, PAR > 0) %>% stats::na.omit()

    # get priors:
    #priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),data = df %>% filter(PAR > 0), family = poisson())

    base::try(model.brms <- brms::brm(brms::bf(equation, ax+alpha+Icomp ~ 1, nl = TRUE), prior = priors.lrc,
                                      data = df.sub, iter = iterations, cores = 3, chains = 1, backend = "cmdstanr"), silent = F)
    # Here it should be df.sub
    base::print(model.brms)

    base::try(model.brms.df <- summary(model.brms)$fixed, silent = T)
    base::try(model.brms.df.ax <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'ax_Intercept'), silent = F)
    base::try(model.brms.df.alpha <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'alpha_Intercept'), silent = F)
    base::try(model.brms.df.Icomp <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'Icomp_Intercept'), silent = F)

    samples <- df.sub %>% dplyr::filter(idx == i) %>% dplyr::select(nee) %>% stats::na.omit() %>% base::nrow()

    baseline <- base::as.Date(base::paste(i, '-01', sep = "")) %>% lubridate::days_in_month()*48 %>% base::as.numeric()

    base::try(results <- base::data.frame(idx = i,
                                          ax.mean = model.brms.df.ax$Estimate,
                                          ax.se = model.brms.df.ax$Est.Error,
                                          ax.Bulk_ESS = model.brms.df.ax$Bulk_ESS,
                                          ax.Tail_ESS = model.brms.df.ax$Tail_ESS,
                                          ax.Rhat = model.brms.df.ax$Rhat,

                                          alpha.mean = model.brms.df.alpha$Estimate,
                                          alpha.se = model.brms.df.alpha$Est.Error,
                                          alpha.Bulk_ESS = model.brms.df.alpha$Bulk_ESS,
                                          alpha.Tail_ESS = model.brms.df.alpha$Tail_ESS,
                                          alpha.Rhat = model.brms.df.alpha$Rhat,

                                          Icomp.mean = model.brms.df.Icomp$Estimate,
                                          Icomp.se = model.brms.df.Icomp$Est.Error,
                                          Icomp.Bulk_ESS = model.brms.df.Icomp$Bulk_ESS,
                                          Icomp.Tail_ESS = model.brms.df.Icomp$Tail_ESS,
                                          Icomp.Rhat = model.brms.df.Icomp$Rhat,
                                          samples= samples/baseline*100), silent = T)

    base::message('YOU DID IT!')

    base::print(results)

    base::try(parms <- parms %>% base::rbind(results), silent = T)

    base::message('Just keep swimming')
  }

  return(parms)
}
# EOF
