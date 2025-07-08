#' @title Fit a Temperature Response Curve (TRC) by an index to get a parameter file
#'
#' @description
#' This function uses the equation:
#' \deqn{\text{NEE} \sim r_0 \cdot \exp \left( \alpha \cdot T + \beta \cdot T^2 \right)}
#'
#' Where NEE is nighttime NEE, that is, ecosystem respiration (\eqn{\mu}mol CO2 m-2 s-1),
#' TA is air temperature (°C),
#' \eqn{r_0} is ecosystem respiration at TA = 0 °C,
#' and alpha and beta are coefficients.
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
#'
#' @param data.frame (dataframe) A dataframe that contains net ecosystem exchange (NEE), an index, air temperature (TA), and photosynthetically active radiation (PAR).
#' @param iterations (numeric) The number of iterations to run `brms::brm()`.
#' @param priors.trc (brmsprior dataframe) The priors for `brms::brm()` to use.
#' Default priors are as follows:
#' ```
#' brms::prior("normal(1, 5)", nlpar = "r0", lb = 0, ub = 5) +
#' brms::prior("normal(0.07, 0.1)", nlpar = "alpha", lb = 0, ub = 0.2) +
#' brms::prior("normal(0, 0.05)", nlpar = "beta", lb = -0.1, ub = 0.1)
#' ```
#' @param idx.colname (character) The name of the column containing the index.
#' @param NEE.colname (character) The name of the column containing NEE.
#' @param PAR.colname (character) The name of the column containing PAR.
#' @param TA.colname (character) The name of the column containing air temperature.
#'
#' @returns (dataframe) Dataframe of parameter values by the index used to fit them.
#' @importFrom magrittr %>%
#' @export
#'
#' @examplesIf !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
#' # Import flux tower data
#' tower.data <- read.csv(system.file("extdata", "AMF_US-Skr_BASE_HH_2-5_Formatted.csv",
#'                                    package = "CarbonExchangeParameters"))
#'
#' # Fit curve parameters for each YearMon:
#' Example_TRC_PARMS_02 <- TRC_PARMS_02(data.frame = tower.data,
#'                                      iterations = 1000,
#'                                      priors.trc = brms::prior("normal(1, 5)",
#'                                                     nlpar = "r0", lb = 0, ub = 5) +
#'                                                   brms::prior("normal(0.07, 0.1)",
#'                                                     nlpar = "alpha", lb = 0, ub = 0.2) +
#'                                                   brms::prior("normal(0, 0.05)",
#'                                                     nlpar = "beta", lb = -0.1, ub = 0.1),
#'                                      idx.colname = 'YearMon',
#'                                      NEE.colname = 'NEE_PI',
#'                                      PAR.colname = 'SW_IN',
#'                                      TA.colname = 'TA_1_1_1')
#'
#'
TRC_PARMS_02 <- function(data.frame = NULL,
                         iterations = NULL,
                         priors.trc = brms::prior("normal(1, 5)", nlpar = "r0", lb = 0, ub = 5) +
                           brms::prior("normal(0.07, 0.1)", nlpar = "alpha", lb = 0, ub = 0.2) +
                           brms::prior("normal(0, 0.05)", nlpar = "beta", lb = -0.1, ub = 0.1),
                         idx.colname = NULL,
                         NEE.colname = NULL,
                         PAR.colname = NULL,
                         TA.colname = NULL){

  # Squelch visible bindings note
  nee <- idx <- TA <- PAR <- NULL

  data.frame$nee <- data.frame[,NEE.colname]
  data.frame$idx <- data.frame[,idx.colname]
  data.frame$PAR <- data.frame[,PAR.colname]
  data.frame$TA <- data.frame[,TA.colname]

  df <- data.frame %>% dplyr::select(idx, nee, PAR, TA)

  if(c("idx") %in% base::names(df)) {
    base::print("GREAT JOB! your dataframe contains idx")
  } else {base::print("The dataframe must include: idx, nee, TA, and PAR")}

  if(c("nee") %in% base::names(df)) {
    base::print("YIPEE! your dataframe contains nee")
  } else {base::print("The dataframe must include: idx, nee, TA, and PAR")}

  if(c("PAR") %in% base::names(df)) {
    base::print("Hooray! your dataframe contains PAR")
  } else {base::print("The dataframe must include: idx, nee, TA, and PAR")}

  if(c("TA") %in% base::names(df)) {
    base::print("Hooray! your dataframe contains TA")
  } else {base::print("The dataframe must include: idx, nee, TA, and PAR")}

  base::try(equation <- nee ~ r0 * exp(alpha*TA + beta*TA*TA), silent = T)

  # PARM Dataframe:
  parms <- base::data.frame(idx = base::as.character(),
                            r0.mean = base::as.numeric(),
                            r0.se = base::as.numeric(),
                            r0.Bulk_ESS = base::as.numeric(),
                            r0.Tail_ESS = base::as.numeric(),
                            r0.Rhat = base::as.numeric(),

                            alpha.mean = base::as.numeric(),
                            alpha.se = base::as.numeric(),
                            alpha.Bulk_ESS = base::as.numeric(),
                            alpha.Tail_ESS = base::as.numeric(),
                            alpha.Rhat = base::as.numeric(),

                            beta.mean = base::as.numeric(),
                            beta.se = base::as.numeric(),
                            beta.Bulk_ESS = base::as.numeric(),
                            beta.Tail_ESS = base::as.numeric(),
                            beta.Rhat = base::as.numeric(),
                            samples = base::as.numeric())

  base::message(" Your dataframe looks good and you are now ready to start fitting models")

  for (i in base::unique(df$idx)){
    base::print(i)

    # Subset the file:
    base::try(df.sub <- df %>% dplyr::filter(idx == i, PAR < 10), silent = T)
    # get priors:

    # priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),
    #                    data = df %>% filter(PAR > 0),
    #                    family = poisson())
    base::try(model.brms <- brms::brm(brms::bf(equation, r0+alpha+beta ~ 1, nl = TRUE),
                                      prior = priors.trc, data = df.sub,
                                      backend = "cmdstanr", iter = iterations, cores = 4, seed = 101), silent = F)
    # Here it should be df.sub
    base::print(model.brms)

    base::try(model.brms.df <- summary(model.brms)$fixed, silent = T)
    base::try(model.brms.df.r0 <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'r0_Intercept'), silent = F)
    base::try(model.brms.df.alpha <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'alpha_Intercept'), silent = F)
    base::try(model.brms.df.beta <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'beta_Intercept'), silent = F)

    samples <- df.sub %>% dplyr::filter(idx == i) %>% dplyr::select(nee) %>% stats::na.omit() %>% base::nrow()

    baseline <- base::as.Date(base::paste(i, '-01', sep = "")) %>% lubridate::days_in_month()*48 %>% base::as.numeric()

    base::try(results <- base::data.frame(idx = i,
                                          r0.mean = model.brms.df.r0$Estimate,
                                          r0.se = model.brms.df.r0$Est.Error,
                                          r0.Bulk_ESS = model.brms.df.r0$Bulk_ESS,
                                          r0.Tail_ESS = model.brms.df.r0$Tail_ESS,
                                          r0.Rhat = model.brms.df.r0$Rhat,

                                          alpha.mean = model.brms.df.alpha$Estimate,
                                          alpha.se = model.brms.df.alpha$Est.Error,
                                          alpha.Bulk_ESS = model.brms.df.alpha$Bulk_ESS,
                                          alpha.Tail_ESS = model.brms.df.alpha$Tail_ESS,
                                          alpha.Rhat = model.brms.df.alpha$Rhat,

                                          beta.mean = model.brms.df.beta$Estimate,
                                          beta.se = model.brms.df.beta$Est.Error,
                                          beta.Bulk_ESS = model.brms.df.beta$Bulk_ESS,
                                          beta.Tail_ESS = model.brms.df.beta$Tail_ESS,
                                          beta.Rhat = model.brms.df.beta$Rhat,
                                          samples = samples/baseline*100), silent = T)

    base::message('YOU DID IT!')
    base::print(results)
    base::try(parms <- parms %>% base::rbind(results), silent = T)

  }

  return(parms)
}
# EOF
