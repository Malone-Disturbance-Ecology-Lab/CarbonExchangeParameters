#' @title Fit a Temperature Response Curve (TRC) by an index to get a parameter file
#'
#' @description
#' This function uses the equation:
#' \deqn{\text{NEE} \sim r_b \cdot \exp \left( E_0 \left( \left( \frac{1}{T_{\text{ref}} - T_0} \right) - \left( \frac{1}{T_{\text{air}} - T_0} \right) \right) \right)}
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
#' @param data.frame (dataframe) A dataframe that contains net ecosystem exchange (NEE), an index, air temperature (TA), and photosynthetically active radiation (PAR).
#' @param iterations (numeric) The number of iterations to run `brms::brm()`.
#' @param priors.trc (brmsprior dataframe) The priors for `brms::brm()` to use.
#' Default priors are as follows:
#' ```
#' brms::prior("normal(0.5, 0.3)", nlpar = "E0", lb = 0.01, ub = 1) +
#' brms::prior("normal(1.0, 0.3)", nlpar = "Rref", lb = 0.01, ub = 1)
#' ```
#' @param idx.colname (character) The name of the column containing the index.
#' @param nee.colname (character) The name of the column containing NEE.
#' @param TA.colname (character) The name of the column containing air temperature.
#' @param PAR.colname (character) The name of the column containing PAR.
#'
#' @returns (dataframe) Dataframe of parameter values by the index used to fit them.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Write examples here
#'
TRC_PARMS_01 <- function(data.frame = NULL,
                         iterations = NULL,
                         priors.trc = brms::prior("normal(0.5, 0.3)", nlpar = "E0", lb = 0.01, ub = 1) +
                           brms::prior("normal(1.0, 0.3)", nlpar = "Rref", lb = 0.01, ub = 1),
                         idx.colname = NULL,
                         nee.colname = NULL,
                         TA.colname = NULL,
                         PAR.colname = NULL){

  # Squelch visible bindings note
  nee <- idx <- TA <- PAR <- NULL

  data.frame$nee <- data.frame[,nee.colname]
  data.frame$idx <- data.frame[,idx.colname]
  data.frame$TA <- data.frame[,TA.colname]
  data.frame$PAR <- data.frame[,PAR.colname]

  df <- data.frame %>% dplyr::select(idx, nee, TA, PAR)

  if(c("idx") %in% base::names(df)) {
    base::print("GREAT JOB! your dataframe contains idx")
  } else {base::print("The dataframe must include: idx, nee, PAR, and TA")}

  if(c("nee") %in% base::names(df)) {
    base::print("YIPEE! your dataframe contains nee")
  } else {base::print("The dataframe must include: idx, nee, PAR, and TA")}

  if(c("PAR") %in% base::names(df)) {
    base::print("Hooray! your dataframe contains PAR")
  } else {base::print("The dataframe must include: idx, nee, PAR, and TA")}

  if(c("TA") %in% base::names(df)) {
    base::print("Hooray! your dataframe contains TA")
  } else {base::print("The dataframe must include: idx, nee, PAR, and TA")}

  base::try(equation <- NEE ~ Rref * exp(E0 * (1 / (10 + 46.02) - 1 / (TA + 46.02))), silent = T)

  # PARM Dataframe:
  base::try(parms <- base::data.frame(idx = base::as.character(),
                                      E0.mean = base::as.numeric(),
                                      E0.se = base::as.numeric(),
                                      E0.Bulk_ESS = base::as.numeric(),
                                      E0.Tail_ESS = base::as.numeric(),
                                      E0.Rhat = base::as.numeric(),

                                      Rref.mean = base::as.numeric(),
                                      Rref.se = base::as.numeric(),
                                      Rref.Bulk_ESS = base::as.numeric(),
                                      Rref.Tail_ESS = base::as.numeric(),
                                      Rref.Rhat = base::as.numeric(),
                                      samples = base::as.numeric()), silent = T)

  base::message("Your dataframe looks good and you are now ready to start fitting models")

  for (i in base::unique(df$idx)){
    base::print(i)

    # Subset the file:
    df.sub <- base::try(df %>% dplyr::filter(idx == i, PAR < 10), silent = TRUE) ## to get night time data

    # get priors:
    base::print(df.sub)
    # priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),
    #                    data = df %>% filter(PAR > 0),
    #                    family = poisson())
    base::try(model.brms <- brms::brm(brms::bf(nee ~ Rref * exp(E0 * (1 / (10 + 46.02) - 1 / (TA + 46.02))), E0+Rref ~ 1, nl = TRUE),
                                      prior = priors.trc, data = df.sub,
                                      backend = "cmdstanr", iter = iterations, cores = 4, seed = 101), silent = F)

    base::print(model.brms)

    base::try(model.brms.df <- summary(model.brms)$fixed, silent = F)

    # Extract parameter estimates correctly
    model.brms.df.E0 <- model.brms.df %>% dplyr::filter(base::rownames(model.brms.df) == "E0_Intercept")
    model.brms.df.Rref <- model.brms.df %>% dplyr::filter(base::rownames(model.brms.df) == "Rref_Intercept")

    base::try(samples <- df.sub %>% dplyr::filter(idx == i) %>% dplyr::select(nee) %>% stats::na.omit() %>% base::nrow(), silent = F)
    base::try(baseline <- base::as.Date(base::paste(i, '-01', sep = "")) %>% lubridate::days_in_month()*48 %>% base::as.numeric(), silent = F)


    base::try(results <- base::data.frame(idx = i,
                                          E0.mean = model.brms.df.E0$Estimate,
                                          E0.se = model.brms.df.E0$Est.Error,
                                          E0.Bulk_ESS = model.brms.df.E0$Bulk_ESS,
                                          E0.Tail_ESS = model.brms.df.E0$Tail_ESS,
                                          E0.Rhat = model.brms.df.E0$Rhat,

                                          Rref.mean = model.brms.df.Rref$Estimate,
                                          Rref.se = model.brms.df.Rref$Est.Error,
                                          Rref.Bulk_ESS = model.brms.df.Rref$Bulk_ESS,
                                          Rref.Tail_ESS = model.brms.df.Rref$Tail_ESS,
                                          Rref.Rhat = model.brms.df.Rref$Rhat,
                                          samples = samples/baseline *100), silent = T)

    base::message('YOU DID IT!')
    base::print(summary(model.brms))
    #print(results)
    base::try(parms <- parms %>% base::rbind(results), silent = T)

  }

  return(parms)
}
#EOF
