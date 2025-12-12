#' @title Fit a Light Response Curve (LRC) by an index to get a parameter file
#'
#' @description
#' This function uses the equation:
#' \deqn{\text{NEE} \sim -a_1 \cdot \exp(-B \cdot \text{PAR}) - Y \cdot \exp(Z \cdot \text{PAR})}
#'
#' Where B and Y are correction factors in which
#' B is the photoinhibitation item,
#' Y is the light saturation item in which Y = a/Pmax,
#' and Z is a correction factor not well defined within the paper.
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
#' brms::prior("normal(-0.12, 0.1)", nlpar = "a", lb = -0.2, ub = 0) +
#' brms::prior("normal(0, 1)", nlpar = "B", lb = -1, ub = 1) +
#' brms::prior("normal(0.2, 0.1)", nlpar = "Y", lb = 0, ub = 0.3) +
#' brms::prior("normal(0.0, 0.1)", nlpar = "Z", lb = -0.1, ub = 0.2)
#' ```
#' @param idx.colname (character) The name of the column containing the index.
#' @param NEE.colname (character) The name of the column containing NEE.
#' @param PAR.colname (character) The name of the column containing PAR.
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
#' Example_LRC_PARMS_05 <- LRC_PARMS_05(data.frame = tower.data,
#'                                      iterations = 1000,
#'                                      priors.lrc = brms::prior("normal(-0.12, 0.1)",
#'                                                     nlpar = "a", lb = -0.2, ub = 0) +
#'                                                   brms::prior("normal(0, 1)",
#'                                                     nlpar = "B", lb = -1, ub = 1) +
#'                                                   brms::prior("normal(0.2, 0.1)",
#'                                                     nlpar = "Y", lb = 0, ub = 0.3) +
#'                                                   brms::prior("normal(0.0, 0.1)",
#'                                                     nlpar = "Z", lb = -0.1, ub = 0.2),
#'                                      idx.colname = 'YearMon',
#'                                      NEE.colname = 'NEE_PI',
#'                                      PAR.colname = 'SW_IN')
#'
#'
LRC_PARMS_05 <- function(data.frame = NULL,
                         iterations = NULL,
                         priors.lrc = brms::prior("normal(-0.12, 0.1)", nlpar = "a", lb = -0.2, ub = 0) +
                           brms::prior("normal(0, 1)", nlpar = "B", lb = -1, ub = 1) +
                           brms::prior("normal(0.2, 0.1)", nlpar = "Y", lb = 0, ub = 0.3) +
                           brms::prior("normal(0.0, 0.1)", nlpar = "Z", lb = -0.1, ub = 0.2),
                         idx.colname = NULL,
                         NEE.colname = NULL,
                         PAR.colname = NULL){

  # Error out if no data frame is provided
  if (base::is.null(data.frame)) stop("Please provide a data frame")

  # Error out if no number of iterations is provided
  if (base::is.null(iterations)) stop("Please provide the number of iterations")

  # Error out if no idx column is provided
  if (base::is.null(idx.colname)) {
    stop("Please provide the name of the column containing the index")
  } else (base::message("GREAT JOB! your dataframe contains an index column"))

  # Error out if no NEE column is provided
  if (base::is.null(NEE.colname)) {
    stop("Please provide the name of the column containing NEE")
  } else (base::message("YIPEE! your dataframe contains NEE"))

  # Error out if no PAR column is provided
  if (base::is.null(PAR.colname)) {
    stop("Please provide the name of the column containing PAR")
  } else (base::message("Hooray! your dataframe contains PAR"))

  # Squelch visible bindings note
  nee <- idx <- PAR <- NULL

  data.frame$nee <- data.frame[,NEE.colname]
  data.frame$idx <- data.frame[,idx.colname]
  data.frame$PAR <- data.frame[,PAR.colname]

  df <- data.frame %>%
    dplyr::select(idx, nee, PAR)

  equation <- nee ~ a * exp(-B * PAR) - Y * exp(-Z * PAR)

  # PARM Dataframe:
  parms <- base::data.frame(idx = base::as.character(),
                            a.mean = base::as.numeric(),
                            a.se = base::as.numeric(),
                            a.Bulk_ESS = base::as.numeric(),
                            a.Tail_ESS = base::as.numeric(),
                            a.Rhat = base::as.numeric(),

                            B.mean = base::as.numeric(),
                            B.se = base::as.numeric(),
                            B.Bulk_ESS = base::as.numeric(),
                            B.Tail_ESS = base::as.numeric(),
                            B.Rhat = base::as.numeric(),

                            Y.mean = base::as.numeric(),
                            Y.se = base::as.numeric(),
                            Y.Bulk_ESS = base::as.numeric(),
                            Y.Tail_ESS = base::as.numeric(),
                            Y.Rhat = base::as.numeric(),

                            Z.mean = base::as.numeric(),
                            Z.se = base::as.numeric(),
                            Z.Bulk_ESS = base::as.numeric(),
                            Z.Tail_ESS = base::as.numeric(),
                            Z.Rhat = base::as.numeric(),
                            samples = base::as.numeric())

  base::message(" Your dataframe looks good and you are now ready to start fitting models")

  for (i in base::unique(df$idx)){
    base::print(i)

    # Subset the file:
    df.sub <- df %>% dplyr::filter(idx == i, PAR > 0) %>% stats::na.omit()

    # get priors:
    #priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),data = df %>% filter(PAR > 0), family = poisson())

    base::try(model.brms <- brms::brm(brms::bf(equation, a+B+Y+Z ~ 1, nl = TRUE),
                                      prior = priors.lrc, data = df.sub, iter = iterations, cores = 3, chains = 1, backend = "cmdstanr"), silent = F)

    base::print(model.brms)

    base::try(model.brms.df <- summary(model.brms)$fixed, silent = T)
    base::try(model.brms.df.a <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'a_Intercept'), silent = F)
    base::try(model.brms.df.B <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'B_Intercept'), silent = F)
    base::try(model.brms.df.Y <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'Y_Intercept'), silent = F)
    base::try(model.brms.df.Z <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'Z_Intercept'), silent = F)

    samples <- df.sub %>% dplyr::filter(idx == i) %>% dplyr::select(nee) %>% stats::na.omit() %>% base::nrow()

    baseline <- base::as.Date(base::paste(i, '-01', sep = "")) %>% lubridate::days_in_month()*48 %>% base::as.numeric()


    base::try(results <- base::data.frame(idx = i,
                                          a.mean = model.brms.df.a$Estimate,
                                          a.se = model.brms.df.a$Est.Error,
                                          a.Bulk_ESS = model.brms.df.a$Bulk_ESS,
                                          a.Tail_ESS = model.brms.df.a$Tail_ESS,
                                          a.Rhat = model.brms.df.a$Rhat,

                                          B.mean = model.brms.df.B$Estimate,
                                          B.se = model.brms.df.B$Est.Error,
                                          B.Bulk_ESS = model.brms.df.B$Bulk_ESS,
                                          B.Tail_ESS = model.brms.df.B$Tail_ESS,
                                          B.Rhat = model.brms.df.B$Rhat,

                                          Y.mean = model.brms.df.Y$Estimate,
                                          Y.se = model.brms.df.Y$Est.Error,
                                          Y.Bulk_ESS = model.brms.df.Y$Bulk_ESS,
                                          Y.Tail_ESS = model.brms.df.Y$Tail_ESS,
                                          Y.Rhat = model.brms.df.Y$Rhat,

                                          Z.mean = model.brms.df.Z$Estimate,
                                          Z.se = model.brms.df.Z$Est.Error,
                                          Z.Bulk_ESS = model.brms.df.Z$Bulk_ESS,
                                          Z.Tail_ESS = model.brms.df.Z$Tail_ESS,
                                          Z.Rhat = model.brms.df.Z$Rhat,
                                          samples= samples/baseline*100), silent = T)

    base::message('YOU DID IT!')
    base::print(results)

    base::try(parms <- parms %>% base::rbind(results), silent = T)
  }

  return(parms)
}
# EOF
