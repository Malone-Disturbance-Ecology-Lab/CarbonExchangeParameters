#' @title Fit a Temperature Response Curve (TRC) by an index to get a parameter file
#'
#' @description
#' This function uses the equation:
#' \deqn{\text{NEE} \sim R_{\text{ref}} \cdot \exp \left( \frac{E_a}{R} \left( \frac{1}{T_{\text{ref}}} - \frac{1}{T} \right) \right) }
#'
#' Where \eqn{R_{ref}} is the reference respiration rate \[\eqn{\mu}mol CO2 m-2 s-1\],
#' \eqn{E_a} is the activation energy \[J/mol\],
#' \eqn{R} is the Universal Gas Constant \[8.314 J mol-1 K-1\],
#' \eqn{T_{ref}} is the reference temperature (25 C or 298K) \[K\],
#' and TA is temperature \[K\].
#'
#' WARNING: TA must be in kelvin.
#'
#' The equation requires air temperature (TA) in kelvin, photosynthetically active radiation (PAR) in \eqn{\mu}mol m-2 s-1, and net ecosystem exchange (NEE) in \eqn{\mu}mol m-2 s-1.
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
#' @param data.frame (data.frame) A data frame that contains net ecosystem exchange (NEE), an index, air temperature (TA), and photosynthetically active radiation (PAR).
#' @param iterations (numeric) The number of iterations to run `brms::brm()`.
#' @param priors.trc (brmsprior data.frame) The priors for `brms::brm()` to use.
#' Default priors are as follows:
#' ```
#' brms::prior("normal(0.5, 0.3)", nlpar = "Ea", lb = 0.01, ub = 1) +
#' brms::prior("normal(.5, .3)", nlpar = "Rref", lb = 0.01, ub = 1)
#' ```
#' @param idx.colname (character) The name of the column containing the index.
#' @param NEE.colname (character) The name of the column containing NEE.
#' @param PAR.colname (character) The name of the column containing PAR.
#' @param TA.colname (character) The name of the column containing air temperature.
#'
#' @returns (data.frame) Data frame of parameter values by the index used to fit them.
#' @importFrom magrittr %>%
#' @export
#'
#' @examplesIf !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
#' # Import flux tower data
#' tower.data <- read.csv(system.file("extdata", "AMF_US-Skr_BASE_HH_2-5_Formatted.csv",
#'                                    package = "CarbonExchangeParameters"))
#'
#' # Fit curve parameters for each YearMon:
#' Example_TRC_PARMS_06 <- TRC_PARMS_06(data.frame = tower.data,
#'                                      iterations = 5000,
#'                                      priors.trc = brms::prior("normal(0.5, 0.3)",
#'                                                     nlpar = "Ea", lb = 0.01, ub = 1) +
#'                                                   brms::prior("normal(.5, .3)",
#'                                                     nlpar = "Rref", lb = 0.01, ub = 1),
#'                                      idx.colname = 'YearMon',
#'                                      NEE.colname = 'NEE_PI',
#'                                      PAR.colname = 'SW_IN',
#'                                      TA.colname = 'TA_1_1_1')
#'
#'
TRC_PARMS_06 <- function(data.frame = NULL,
                         iterations = NULL,
                         priors.trc = brms::prior("normal(0.5, 0.3)", nlpar = "Ea", lb = 0.01, ub = 1) +
                           brms::prior("normal(.5, .3)", nlpar = "Rref", lb = 0.01, ub = 1),
                         idx.colname = NULL,
                         NEE.colname = NULL,
                         PAR.colname = NULL,
                         TA.colname = NULL){

  # Error out if no data frame is provided
  if (base::is.null(data.frame)) stop("Please provide a data frame")

  # Error out if no number of iterations is provided
  if (base::is.null(iterations)) stop("Please provide the number of iterations")

  # Error out if no idx column is provided
  if (base::is.null(idx.colname)) {
    stop("Please provide the name of the column containing the index")
  } else (base::message("GREAT JOB! your data frame contains an index column"))

  # Error out if no NEE column is provided
  if (base::is.null(NEE.colname)) {
    stop("Please provide the name of the column containing NEE")
  } else (base::message("YIPEE! your data frame contains NEE"))

  # Error out if no TA column is provided
  if (base::is.null(TA.colname)) {
    stop("Please provide the name of the column containing TA")
  } else (base::message("Yay! your data frame contains TA"))

  # Error out if no PAR column is provided
  if (base::is.null(PAR.colname)) {
    stop("Please provide the name of the column containing PAR")
  } else (base::message("Hooray! your data frame contains PAR"))

  # Squelch visible bindings note
  nee <- idx <- TA <- PAR <- NULL

  data.frame$nee <- data.frame[,NEE.colname]
  data.frame$idx <- data.frame[,idx.colname]
  data.frame$PAR <- data.frame[,PAR.colname]
  data.frame$TA <- data.frame[,TA.colname]

  df <- data.frame %>% dplyr::select(idx, nee, PAR, TA)

  base::try(equation <- nee ~ Rref * exp((Ea/R)*((1/Tref)-(1/TA))), silent = T)

  # PARM Data frame:
  base::try(parms <- base::data.frame(idx = base::as.character(),
                                      Ea.mean = base::as.numeric(),
                                      Ea.se = base::as.numeric(),
                                      Ea.Bulk_ESS = base::as.numeric(),
                                      Ea.Tail_ESS = base::as.numeric(),
                                      Ea.Rhat = base::as.numeric(),

                                      Rref.mean = base::as.numeric(),
                                      Rref.se = base::as.numeric(),
                                      Rref.Bulk_ESS = base::as.numeric(),
                                      Rref.Tail_ESS = base::as.numeric(),
                                      Rref.Rhat = base::as.numeric(),
                                      samples = base::as.numeric()), silent = T)

  base::message(" Your data frame looks good and you are now ready to start fitting models")

  for (i in base::unique(df$idx)){
    base::print(i)

    # Subset the file:
    df.sub <- base::try(df %>% dplyr::filter(idx == i, PAR < 10), silent = TRUE)

    # get priors:
    base::print(df.sub)
    # priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),
    #                    data = df %>% filter(PAR > 0),
    #                    family = poisson())
    base::try(model.brms <- brms::brm(brms::bf(nee ~ Rref * exp((Ea/8.314)*((1.0/298)-(1/TA))), Ea+Rref ~ 1, nl = TRUE),
                                      prior = priors.trc, data = df.sub,
                                      backend = "cmdstanr", iter = iterations, cores = 4, seed = 101), silent = F)

    base::print(model.brms)

    base::try(model.brms.df <- summary(model.brms)$fixed, silent = F)

    base::try(model.brms.df.Ea <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'Ea_Intercept'), silent = F)
    base::try(model.brms.df.Rref <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'Rref_Intercept'), silent = F)

    base::try(samples <- df.sub %>% dplyr::filter(idx == i) %>% dplyr::select(nee) %>% stats::na.omit() %>% base::nrow(), silent = F)
    base::try(baseline <- base::as.Date(base::paste(i, '-01', sep = "")) %>% lubridate::days_in_month()*48 %>% base::as.numeric(), silent = F)


    base::try(results <- base::data.frame(idx = i,
                                          Ea.mean = model.brms.df.Ea$Estimate,
                                          Ea.se = model.brms.df.Ea$Est.Error,
                                          Ea.Bulk_ESS = model.brms.df.Ea$Bulk_ESS,
                                          Ea.Tail_ESS = model.brms.df.Ea$Tail_ESS,
                                          Ea.Rhat = model.brms.df.Ea$Rhat,

                                          Rref.mean = model.brms.df.Rref$Estimate,
                                          Rref.se = model.brms.df.Rref$Est.Error,
                                          Rref.Bulk_ESS = model.brms.df.Rref$Bulk_ESS,
                                          Rref.Tail_ESS = model.brms.df.Rref$Tail_ESS,
                                          Rref.Rhat = model.brms.df.Rref$Rhat,
                                          samples= samples/baseline*100), silent = T)

    base::message('YOU DID IT!')
    base::print(results)
    base::try(parms <- parms %>% base::rbind(results), silent = T)

  }

  return(parms)
}
#EOF
