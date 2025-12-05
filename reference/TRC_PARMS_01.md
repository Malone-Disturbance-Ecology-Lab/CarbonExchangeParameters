# Fit a Temperature Response Curve (TRC) by an index to get a parameter file

This function uses the equation: \$\$\text{NEE} \sim r_b \cdot \exp
\left( E_0 \left( \left( \frac{1}{T\_{\text{ref}} - T_0} \right) -
\left( \frac{1}{T\_{\text{air}} - T_0} \right) \right) \right)\$\$

The equation requires air temperature (TA) in degrees Celsius,
photosynthetically active radiation (PAR) in \\\mu\\mol m-2 s-1, and net
ecosystem exchange (NEE) in \\\mu\\mol m-2 s-1.

## Usage

``` r
TRC_PARMS_01(
  data.frame = NULL,
  iterations = NULL,
  priors.trc = brms::prior("normal(0.5, 0.3)", nlpar = "E0", lb = 0.01, ub = 1) +
    brms::prior("normal(1.0, 0.3)", nlpar = "Rref", lb = 0.01, ub = 1),
  idx.colname = NULL,
  NEE.colname = NULL,
  TA.colname = NULL,
  PAR.colname = NULL
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

      brms::prior("normal(0.5, 0.3)", nlpar = "E0", lb = 0.01, ub = 1) +
      brms::prior("normal(1.0, 0.3)", nlpar = "Rref", lb = 0.01, ub = 1)

- idx.colname:

  (character) The name of the column containing the index.

- NEE.colname:

  (character) The name of the column containing NEE.

- TA.colname:

  (character) The name of the column containing air temperature.

- PAR.colname:

  (character) The name of the column containing PAR.

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
Example_TRC_PARMS_01 <- TRC_PARMS_01(data.frame = tower.data,
                                     iterations = 1000,
                                     priors.trc = brms::prior("normal(0.5, 0.3)",
                                                    nlpar = "E0", lb = 0.01, ub = 1) +
                                                  brms::prior("normal(1.0, 0.3)",
                                                    nlpar = "Rref", lb = 0.01, ub = 1),
                                     idx.colname = 'YearMon',
                                     NEE.colname = 'NEE_PI',
                                     PAR.colname = 'SW_IN',
                                     TA.colname = 'TA_1_1_1')

}
```
