#' @title Fit a Temperature Response Curve (TRC) by an index to get a parameter file
#'
#' @description
#' This function uses the equation:
#' \deqn{\text{NEE} \sim a * \exp \left(b*T\right)}
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
#' brms::prior(normal(0.2, 1), nlpar = "a", lb = 0.1, ub = 1) +
#' brms::prior(normal(0.5, 0.03), nlpar = "b", lb = 0.001, ub = 0.9)
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
TRC_PARMS_05 <- function(data.frame = NULL,
                         iterations = NULL,
                         priors.trc = brms::prior(normal(0.2 , 1), nlpar = "a", lb = 0.1, ub = 1) +
                           brms::prior(normal(0.5, 0.03), nlpar = "b", lb = 0.001, ub = 0.9),
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
  base::try(equation <- nee ~ a * exp(b*TA), silent = T)

  # PARM Dataframe:
  base::try(parms <- base::data.frame(idx = base::as.character(),
                                      a.mean = base::as.numeric(),
                                      a.se = base::as.numeric(),
                                      a.Bulk_ESS = base::as.numeric(),
                                      a.Tail_ESS = base::as.numeric(),
                                      a.Rhat = base::as.numeric(),

                                      b.mean = base::as.numeric(),
                                      b.se = base::as.numeric(),
                                      b.Bulk_ESS = base::as.numeric(),
                                      b.Tail_ESS = base::as.numeric(),
                                      b.Rhat = base::as.numeric(),
                                      samples = base::as.numeric()), silent = T)

  base::message(" Your dataframe looks good and you are now ready to start fitting models")

  for (i in base::unique(df$idx)){
    base::print(i)

    # Subset the file:
    base::try(df.sub <- df %>% dplyr::filter(idx == i, PAR < 10), silent = T)
    # get priors:

    base::try(model.brms <- brms::brm(brms::bf(nee ~ a * exp(b*TA), a+b ~ 1, nl = TRUE),
                                      prior = priors.trc, data = df.sub,
                                      backend = "cmdstanr", iter = iterations, cores = 4, seed = 101), silent = F)

    base::print(model.brms)

    base::try(model.brms.df <- summary(model.brms)$fixed, silent = F)

    base::try(model.brms.df.a <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'a_Intercept'), silent = F)
    base::try(model.brms.df.b <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'b_Intercept'), silent = F)

    base::try(samples <- df.sub %>% dplyr::filter(idx == i) %>% dplyr::select(nee) %>% stats::na.omit() %>% base::nrow(), silent = F)
    base::try(baseline <- base::as.Date(base::paste(i, '-01', sep = "")) %>% lubridate::days_in_month()*48 %>% base::as.numeric(), silent = F)


    base::try(results <- base::data.frame(idx = i,
                                          a.mean = model.brms.df.a$Estimate,
                                          a.se = model.brms.df.a$Est.Error,
                                          a.Bulk_ESS = model.brms.df.a$Bulk_ESS,
                                          a.Tail_ESS = model.brms.df.a$Tail_ESS,
                                          a.Rhat = model.brms.df.a$Rhat,

                                          b.mean = model.brms.df.b$Estimate,
                                          b.se = model.brms.df.b$Est.Error,
                                          b.Bulk_ESS = model.brms.df.b$Bulk_ESS,
                                          b.Tail_ESS = model.brms.df.b$Tail_ESS,
                                          b.Rhat = model.brms.df.b$Rhat,
                                          samples= samples/baseline*100), silent = T)

    base::message('YOU DID IT!')
    base::print(results)
    base::try(parms <- parms %>% base::rbind(results), silent = T)

  }

  base::return(parms)
}
#EOF
