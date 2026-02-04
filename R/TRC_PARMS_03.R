#' @title Fit a Temperature Response Curve (TRC) by an index to get a parameter file
#'
#' @description  Author Ammara Talib
#' This function uses the equation:
#' \deqn{\text{NEE} \sim a \left( T - T_{\text{opt}} \right)^2 + E_{\text{Rmax}}}
#'
#' Where
#' \eqn{a} = a parameter estimate of the fitted quadratic function,
#' \eqn{T} = temperature in Celsius,
#' \eqn{T_{opt}} = optimum temperature, a parameter, vertex of the parabola,
#' and \eqn{E_{Rmax}} = maximum respiration, a parameter, vertex of the parabola.
#'
#' Reference: https://www.nature.com/articles/s41559-023-02121-w#Equ1
#'
#' The equation requires air temperature (TA) in degrees Celsius, photosynthetically active radiation (PAR) in \eqn{\mu}mol m-2 s-1, and net ecosystem exchange (NEE) in \eqn{\mu}mol m-2 s-1.
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
#' brms::prior("normal(0.08, 0.03)", nlpar = "a", lb = 0, ub = 0.9) +
#' brms::prior("normal(24, 4)", nlpar = "Topt", lb = 10, ub = 40) +
#' brms::prior("normal(3.8, 5)", nlpar = "ERmax", lb = 0, ub = 10)
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
#' Example_TRC_PARMS_03 <- TRC_PARMS_03(data.frame = tower.data,
#'                                      iterations = 5000,
#'                                      priors.trc = brms::prior("normal(0.08, 0.03)",
#'                                                     nlpar = "a", lb = 0, ub = 0.9) +
#'                                                   brms::prior("normal(24, 4)",
#'                                                     nlpar = "Topt", lb = 10, ub = 40) +
#'                                                   brms::prior("normal(3.8, 5)",
#'                                                     nlpar = "ERmax", lb = 0, ub = 10),
#'                                      idx.colname = 'YearMon',
#'                                      NEE.colname = 'NEE_PI',
#'                                      PAR.colname = 'SW_IN',
#'                                      TA.colname = 'TA_1_1_1')
#'
#'
TRC_PARMS_03 <- function(data.frame = NULL,
                         iterations = NULL,
                         priors.trc = brms::prior("normal(0.08, 0.03)", nlpar = "a", lb = 0, ub = 0.9) +
                           brms::prior("normal(24, 4)", nlpar = "Topt", lb = 10, ub = 40) +
                           brms::prior("normal(3.8, 5)", nlpar = "ERmax", lb = 0, ub = 10),
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

  base::try(equation <- nee ~ a * ((TA - Topt)^2) + ERmax, silent = T)

  # PARM Data frame:
  base::try(parms <- base::data.frame(idx = base::as.character(),
                                      a.mean = base::as.numeric(),
                                      a.se = base::as.numeric(),
                                      a.Bulk_ESS = base::as.numeric(),
                                      a.Tail_ESS = base::as.numeric(),
                                      a.Rhat = base::as.numeric(),

                                      Topt.mean = base::as.numeric(),
                                      Topt.se = base::as.numeric(),
                                      Topt.Bulk_ESS = base::as.numeric(),
                                      Topt.Tail_ESS = base::as.numeric(),
                                      Topt.Rhat = base::as.numeric(),

                                      ERmax.mean = base::as.numeric(),
                                      ERmax.se = base::as.numeric(),
                                      ERmax.Bulk_ESS = base::as.numeric(),
                                      ERmax.Tail_ESS = base::as.numeric(),
                                      ERmax.Rhat = base::as.numeric(),
                                      samples = base::as.numeric()), silent = T)

  base::message(" Your data frame looks good and you are now ready to start fitting models")

  for (i in base::unique(df$idx)){
    base::print(i)

    # Subset the file:
    df.sub <- base::try(df %>% dplyr::filter(idx == i, PAR < 10), silent = TRUE)

    # get priors:
    base::print(df.sub)

    base::try(model.brms <- brms::brm(brms::bf(nee ~ a * ((TA - Topt)^2) + ERmax, a+Topt+ERmax ~ 1, nl = TRUE),
                                      prior = priors.trc, data = df.sub,
                                      backend = "cmdstanr", iter = iterations, cores = 4, seed = 101), silent = F)

    base::print(model.brms)

    base::try(model.brms.df <- summary(model.brms)$fixed, silent = F)

    base::try(model.brms.df.a <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'a_Intercept'), silent = F)
    base::try(model.brms.df.Topt <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'Topt_Intercept'), silent = F)
    base::try(model.brms.df.ERmax <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'ERmax_Intercept'), silent = F)

    base::try(samples <- df.sub %>% dplyr::filter(idx == i) %>% dplyr::select(nee) %>% stats::na.omit() %>% base::nrow(), silent = F)
    base::try(baseline <- base::as.Date(base::paste(i, '-01', sep = "")) %>% lubridate::days_in_month()*48 %>% base::as.numeric(), silent = F)


    base::try(results <- base::data.frame(idx = i,
                                          a.mean = model.brms.df.a$Estimate,
                                          a.se = model.brms.df.a$Est.Error,
                                          a.Bulk_ESS = model.brms.df.a$Bulk_ESS,
                                          a.Tail_ESS = model.brms.df.a$Tail_ESS,
                                          a.Rhat = model.brms.df.a$Rhat,

                                          Topt.mean = model.brms.df.Topt$Estimate,
                                          Topt.se = model.brms.df.Topt$Est.Error,
                                          Topt.Bulk_ESS = model.brms.df.Topt$Bulk_ESS,
                                          Topt.Tail_ESS = model.brms.df.Topt$Tail_ESS,
                                          Topt.Rhat = model.brms.df.Topt$Rhat,

                                          ERmax.mean = model.brms.df.ERmax$Estimate,
                                          ERmax.se = model.brms.df.ERmax$Est.Error,
                                          ERmax.Bulk_ESS = model.brms.df.ERmax$Bulk_ESS,
                                          ERmax.Tail_ESS = model.brms.df.ERmax$Tail_ESS,
                                          ERmax.Rhat = model.brms.df.ERmax$Rhat,
                                          samples= samples/baseline*100), silent = T)

    base::message('YOU DID IT!')
    base::print(results)
    try(parms <- parms %>% base::rbind(results), silent = T)

  }

  return(parms)
}
#EOF
