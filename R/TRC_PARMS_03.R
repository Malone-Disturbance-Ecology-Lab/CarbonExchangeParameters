#' @title Fit a Temperature Response Curve (TRC) by an index to get a parameter file
#'
#' @description
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
#' The equation requires air temperature (TA) in degrees Celsius, photosynthetically active radiation (PAR) in 𝜇mol m-2 s-1, and net ecosystem exchange (NEE) in 𝜇mol m-2 s-1.
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
#' brms::prior(normal(0.08, 0.03), nlpar = "a", lb = 0, ub = 0.9) +
#' brms::prior(normal(24, 4), nlpar = "Topt", lb = 10, ub = 40) +
#' brms::prior(normal(3.8, 5), nlpar = "ERmax", lb = 0, ub = 10)
#' ```
#' @param idx.colname (character) The name of the column containing the index.
#' @param NEE.colname (character) The name of the column containing NEE.
#' @param PAR.colname (character) The name of the column containing PAR.
#' @param TA.colname (character) The name of the column containing air temperature.
#'
#' @returns (dataframe) Dataframe of parameter values by the index used to fit them.
#' @export
#'
#' @examples
#' # Write examples here
#'
TRC_PARMS_03 <- function(data.frame = NULL,
                         iterations = NULL,
                         priors.trc = brms::prior(normal(0.08, 0.03), nlpar = "a", lb = 0, ub = 0.9) +
                           brms::prior(normal(24, 4), nlpar = "Topt", lb = 10, ub = 40) +
                           brms::prior(normal(3.8, 5), nlpar = "ERmax", lb = 0, ub = 10),
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
    print("Hooray! your dataframe contains TA, this is a test")
  } else {base::print("The dataframe must include: idx, nee, TA, and PAR")}
  base::try(equation <- nee ~ a * ((TA - Topt)^2) + ERmax, silent = T)


  # PARM Dataframe:
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

  base::message(" Your dataframe looks good and you are now ready to start fitting models")

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

    base::try(samples <- data.frame %>% dplyr::filter(YearMon == i) %>% dplyr::select(nee) %>% stats::na.omit() %>% base::nrow(), silent = F)
    base::try(baseline <- base::as.Date(base::paste(i, '-01', sep = "")) %>% lubridate::days_in_month()*48 %>% base::as.numeric, silent = F)


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
