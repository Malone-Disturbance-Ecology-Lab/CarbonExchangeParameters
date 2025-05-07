#' @title Fit a Light Response Curve (LRC) by an index to get a parameter file
#'
#' @description
#' This function uses the equation:
#' \deqn{\text{NEE} \sim a_1 \left(\frac{1-b \cdot \text{PAR}}{1 + \frac{a_1}{a_x} \cdot \text{PAR}} \right) \cdot \text{PAR} - r}
#'
#' Where \eqn{r} is ecosystem respiration (𝜇mol CO2 m-2 s-1),
#' \eqn{a_1} is the apparent quantum efficiency of CO2 uptake (CO2),
#' \eqn{b} is the photoinhibitation correction factor,
#' and \eqn{a_x} is the maximum CO2 uptake rate on the ecosystem scale.
#'
#' The equation requires photosynthetically active radiation, PAR, in 𝜇mol m-2 s-1 and net ecosystem exchange, NEE, in 𝜇mol m-2 s-1.
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
#' brms::prior(normal(-0.01, 0.1), nlpar = "a1", lb = -0.2, ub = 0) +
#' brms::prior(normal(-7.65, 0.33), nlpar = "ax", lb = -30, ub = -5) +
#' brms::prior(normal(2.10, 0.11), nlpar = "r", lb = 1.9, ub = 2.2) +
#' brms::prior(normal(.3, .2), nlpar = "b", lb = 0, ub= 3)
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
#' # Write examples here
#'
LRC_PARMS_04 <- function(data.frame = NULL,
                         iterations = NULL,
                         priors.lrc = brms::prior(normal(-0.01, 0.1), nlpar = "a1", lb = -0.2, ub = 0) +
                           brms::prior(normal(-7.65, 0.33), nlpar = "ax", lb = -30, ub = -5) +
                           brms::prior(normal(2.10, 0.11), nlpar = "r", lb = 1.9, ub = 2.2) +
                           brms::prior(normal(.3, .2), nlpar = "b", lb = 0, ub = 3),
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

  equation <- nee ~ a1 * ((1-b*PAR) / (1+(a1/ax)*PAR))*PAR - r



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
                            a1.mean = base::as.numeric(),
                            a1.se = base::as.numeric(),
                            a1.Bulk_ESS = base::as.numeric(),
                            a1.Tail_ESS = base::as.numeric(),
                            a1.Rhat = base::as.numeric(),

                            ax.mean = base::as.numeric(),
                            ax.se = base::as.numeric(),
                            ax.Bulk_ESS = base::as.numeric(),
                            ax.Tail_ESS = base::as.numeric(),
                            ax.Rhat = base::as.numeric(),

                            r.mean = base::as.numeric(),
                            r.se = base::as.numeric(),
                            r.Bulk_ESS = base::as.numeric(),
                            r.Tail_ESS = base::as.numeric(),
                            r.Rhat = base::as.numeric(),
                            samples = base::as.numeric(),

                            b.mean = base::as.numeric(),
                            b.se = base::as.numeric(),
                            b.Bulk_ESS = base::as.numeric(),
                            b.Tail_ESS = base::as.numeric(),
                            b.Rhat = base::as.numeric()
  )


  base::message(" Your dataframe looks good and you are now ready to start fitting models")

  for (i in base::unique(df$idx)){
    base::print(i)

    # Subset the file:
    df.sub <- df %>% dplyr::filter(idx == i, PAR > 0) %>% stats::na.omit()

    # get priors:
    #priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),data = df %>% filter(PAR > 0), family = poisson())

    base::try(model.brms <- brms::brm(brms::bf(equation, a1+ax+r+b ~ 1, nl = TRUE),
                                      prior = priors.lrc, data = df, iter = iterations, cores = 3, chains = 1, backend = "cmdstanr"), silent = F)

    base::print(model.brms)

    base::try(model.brms.df <- summary(model.brms)$fixed, silent = T)
    base::try(model.brms.df.a1 <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'a1_Intercept'), silent = F)
    base::try(model.brms.df.ax <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'ax_Intercept'), silent = F)
    base::try(model.brms.df.r <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'r_Intercept'), silent = F)
    base::try(model.brms.df.b <- model.brms.df %>% dplyr::filter(base::row.names(model.brms.df) == 'b_Intercept'), silent = F)

    samples <- df.sub %>% dplyr::filter(idx == i) %>% dplyr::select(nee) %>% stats::na.omit() %>% base::nrow()

    baseline <- base::as.Date(base::paste(i, '-01', sep = "")) %>% lubridate::days_in_month()*48 %>% base::as.numeric()


    base::try(results <- base::data.frame(idx = i,
                                          a1.mean = model.brms.df.a1$Estimate,
                                          a1.se = model.brms.df.a1$Est.Error,
                                          a1.Bulk_ESS = model.brms.df.a1$Bulk_ESS,
                                          a1.Tail_ESS = model.brms.df.a1$Tail_ESS,
                                          a1.Rhat = model.brms.df.a1$Rhat,

                                          ax.mean = model.brms.df.ax$Estimate,
                                          ax.se = model.brms.df.ax$Est.Error,
                                          ax.Bulk_ESS = model.brms.df.ax$Bulk_ESS,
                                          ax.Tail_ESS = model.brms.df.ax$Tail_ESS,
                                          ax.Rhat = model.brms.df.ax$Rhat,

                                          r.mean = model.brms.df.r$Estimate,
                                          r.se = model.brms.df.r$Est.Error,
                                          r.Bulk_ESS = model.brms.df.r$Bulk_ESS,
                                          r.Tail_ESS = model.brms.df.r$Tail_ESS,
                                          r.Rhat = model.brms.df.r$Rhat,

                                          b.mean = model.brms.df.r$Estimate,
                                          b.se = model.brms.df.r$Est.Error,
                                          b.Bulk_ESS = model.brms.df.r$Bulk_ESS,
                                          b.Tail_ESS = model.brms.df.r$Tail_ESS,
                                          b.Rhat = model.brms.df.r$Rhat,

                                          samples = samples/baseline*100), silent = T)

    base::message('YOU DID IT!')

    base::print(results)

    base::try(parms <- parms %>% base::rbind(results), silent = T)

    base::message('Just keep swimming')
  }

  return(parms)
}
# EOF
