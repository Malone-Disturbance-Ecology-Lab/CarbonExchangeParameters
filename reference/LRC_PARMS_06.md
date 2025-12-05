# Fit a Light Response Curve (LRC) by an index to get a parameter file

This function uses the equation: \$\$\text{NEE} \sim a_x \left( 1 - \exp
\left( a_1 \cdot \left( \text{PAR} - I\_{\text{comp}} \right) \right)
\right)\$\$

Where NEE is daytime net ecosystem exchange (\\\mu\\mol CO2 m-2 s-1),
\\a_x\\ is the maximum photosynthesis rate (\\\mu\\mol CO2 m-2 s-1),
\\a_1\\ is a coefficient, and \\I\_{\text{comp}}\\ is the compensation
light level with the same unit of PAR.

The equation requires photosynthetically active radiation, PAR, in
\\\mu\\mol m-2 s-1 and net ecosystem exchange, NEE, in \\\mu\\mol m-2
s-1.

## Usage

``` r
LRC_PARMS_06(
  data.frame = NULL,
  iterations = NULL,
  priors.lrc = brms::prior("normal(-5, 5)", nlpar = "ax", lb = -30, ub = 0) +
    brms::prior("normal(0.1, 1)", nlpar = "alpha", lb = 0, ub = 1) +
    brms::prior("normal(300, 100)", nlpar = "Icomp", lb = 0, ub = 500),
  idx.colname = NULL,
  NEE.colname = NULL,
  PAR.colname = NULL
)
```

## Arguments

- data.frame:

  (dataframe) A dataframe that contains net ecosystem exchange (NEE), an
  index, and photosynthetically active radiation (PAR).

- iterations:

  (numeric) The number of iterations to run
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

- priors.lrc:

  (brmsprior dataframe) The priors for
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) to
  use. Default priors are as follows:

      brms::prior("normal(-5, 5)", nlpar = "ax", lb = -30, ub = 0) +
      brms::prior("normal(0.1, 1)", nlpar = "alpha", lb = 0, ub = 1) +
      brms::prior("normal(300, 100)", nlpar = "Icomp", lb = 0, ub = 500)

- idx.colname:

  (character) The name of the column containing the index.

- NEE.colname:

  (character) The name of the column containing NEE.

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
Example_LRC_PARMS_06 <- LRC_PARMS_06(data.frame = tower.data,
                                     iterations = 1000,
                                     priors.lrc = brms::prior("normal(-5, 5)",
                                                    nlpar = "ax", lb = -30, ub = 0) +
                                                  brms::prior("normal(0.1, 1)",
                                                    nlpar = "alpha", lb = 0, ub = 1) +
                                                  brms::prior("normal(300, 100)",
                                                    nlpar = "Icomp", lb = 0, ub = 500),
                                     idx.colname = 'YearMon',
                                     NEE.colname = 'NEE_PI',
                                     PAR.colname = 'SW_IN')

}
```
