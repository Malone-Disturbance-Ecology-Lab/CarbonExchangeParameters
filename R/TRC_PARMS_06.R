#' @title Fit a Temperature Response Curve (TRC) by an index to get a parameter file
#'
#' @description
#' This function uses the equation:
#' \deqn{\text{NEE} \sim R_{\text{ref}} \cdot \exp \left( \frac{E_a}{R} \left( \frac{1}{T_{\text{ref}}} - \frac{1}{T} \right) \right) }
#'
#' Where \eqn{R_{ref}} is the reference respiration rate \[𝜇mol CO2 m-2 s-1\],
#' \eqn{E_a} is the activation energy \[J/mol\],
#' \eqn{R} is the Universal Gas Constant \[8.314 J mol-1 K-1\],
#' \eqn{T_{ref}} is the reference temperature (25 C or 298K) \[K\],
#' and TA is temperature \[K\].
#'
#' WARNING: TA must be in kelvin.
#'
#' The equation requires air temperature (TA) in kelvin, photosynthetically active radiation (PAR) in 𝜇mol m-2 s-1, and net ecosystem exchange (NEE) in 𝜇mol m-2 s-1.
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
#' @param data.frame (dataframe) A dataframe that contains net ecosystem exchange (NEE), an index, air temperature (TA), and photosynthetically active radiation (PAR).
#' @param iterations (numeric) The number of iterations to run `brms::brm()`.
#' @param priors.trc (brmsprior dataframe) The priors for `brms::brm()` to use.
#' Default priors are as follows:
#' ```
#' brms::prior(normal(0.5, 0.3), nlpar = "Ea", lb = 0.01, ub = 1) +
#' brms::prior(normal(.5, .3), nlpar = "Rref", lb = 0.01, ub = 1)
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
#' @examples
#' # Write examples here
#'
TRC_PARMS_06 <- function(data.frame = NULL,
                         iterations = NULL,
                         priors.trc = brms::prior(normal(0.5, 0.3), nlpar = "Ea", lb = 0.01, ub = 1) +
                           brms::prior(normal(.5, .3), nlpar = "Rref", lb = 0.01, ub = 1),
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
    base::print("Hooray! your dataframe contains TA, this is a test")
  } else {base::print("The dataframe must include: idx, nee, TA, and PAR")}
  base::try(equation <- nee ~ Rref * exp((Ea/R)*((1/Tref)-(1/TA))), silent = T)


  # PARM Dataframe:
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

  base::message(" Your dataframe looks good and you are now ready to start fitting models")

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
