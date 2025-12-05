# Fit a Temperature Response Curve (TRC) by an index to get a parameter file

This function uses the equation: \$\$\text{NEE} \sim a \left( T -
T\_{\text{opt}} \right)^2 + E\_{\text{Rmax}}\$\$

Where \\a\\ = a parameter estimate of the fitted quadratic function,
\\T\\ = temperature in Celsius, \\T\_{opt}\\ = optimum temperature, a
parameter, vertex of the parabola, and \\E\_{Rmax}\\ = maximum
respiration, a parameter, vertex of the parabola.

Reference: https://www.nature.com/articles/s41559-023-02121-w#Equ1

The equation requires air temperature (TA) in degrees Celsius,
photosynthetically active radiation (PAR) in \\\mu\\mol m-2 s-1, and net
ecosystem exchange (NEE) in \\\mu\\mol m-2 s-1.

## Usage

``` r
TRC_PARMS_03(
  data.frame = NULL,
  iterations = NULL,
  priors.trc = brms::prior("normal(0.08, 0.03)", nlpar = "a", lb = 0, ub = 0.9) +
    brms::prior("normal(24, 4)", nlpar = "Topt", lb = 10, ub = 40) +
    brms::prior("normal(3.8, 5)", nlpar = "ERmax", lb = 0, ub = 10),
  idx.colname = NULL,
  NEE.colname = NULL,
  PAR.colname = NULL,
  TA.colname = NULL
)
```

## Arguments

- data.frame:

  (dataframe) A dataframe that contains net ecosystem exchange (NEE), an
  index, air temperature (TA), and photosynthetically active radiation
  (PAR).

- iterations:

  (numeric) The number of iterations to run
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

- priors.trc:

  (brmsprior dataframe) The priors for
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) to
  use. Default priors are as follows:

      brms::prior("normal(0.08, 0.03)", nlpar = "a", lb = 0, ub = 0.9) +
      brms::prior("normal(24, 4)", nlpar = "Topt", lb = 10, ub = 40) +
      brms::prior("normal(3.8, 5)", nlpar = "ERmax", lb = 0, ub = 10)

- idx.colname:

  (character) The name of the column containing the index.

- NEE.colname:

  (character) The name of the column containing NEE.

- PAR.colname:

  (character) The name of the column containing PAR.

- TA.colname:

  (character) The name of the column containing air temperature.

## Value

(dataframe) Dataframe of parameter values by the index used to fit them.

## Details

Model parameters are fit using the R package `brms`.

Rhat (Potential Scale Reduction Factor): Indicates how well the
different Markov chains in your analysis have converged to the same
posterior distribution. Ideally, Rhat should be close to 1 for all
parameters. A high Rhat value suggests potential convergence issues and
the need to run the chains longer.

Bulk ESS (Effective Sample Size - Bulk): Estimates the effective number
of independent samples from the central part of the posterior
distribution.

Tail ESS (Effective Sample Size - Tail): Estimates the effective number
of independent samples from the tails of the posterior distribution.
Important for assessing the reliability of quantile estimates (e.g., 95%
confidence intervals).

Key points to remember: Aim for Rhat close to 1 and high values for both
Bulk ESS and Tail ESS.

## Examples

``` r
if (FALSE) { # !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
# Import flux tower data
tower.data <- read.csv(system.file("extdata", "AMF_US-Skr_BASE_HH_2-5_Formatted.csv",
                                   package = "CarbonExchangeParameters"))

# Fit curve parameters for each YearMon:
Example_TRC_PARMS_03 <- TRC_PARMS_03(data.frame = tower.data,
                                     iterations = 5000,
                                     priors.trc = brms::prior("normal(0.08, 0.03)",
                                                    nlpar = "a", lb = 0, ub = 0.9) +
                                                  brms::prior("normal(24, 4)",
                                                    nlpar = "Topt", lb = 10, ub = 40) +
                                                  brms::prior("normal(3.8, 5)",
                                                    nlpar = "ERmax", lb = 0, ub = 10),
                                     idx.colname = 'YearMon',
                                     NEE.colname = 'NEE_PI',
                                     PAR.colname = 'SW_IN',
                                     TA.colname = 'TA_1_1_1')

}
```
