## Carbon Exchange Parameters

Flux partitioning methods use a light response curve (LRC). We can use daytime NEE (NEEday) from tower sites to fit LRCs and nighttime NEE (NEEnight) to estimate ecosystem respiration with temperature response curves (TRC). Using daytime data,  fit LRCs:

NEE ~ (a1 * PAR * ax)/(a1 * PAR + ax) + r									Eq. 1

Where r is respiration (𝜇mol CO2 m-2 s-1), a1 is the apparent quantum efficiency of CO2 uptake (Reichstein et al., 2012)), PAR is solar radiation (𝜇mol m-2 s-1). ax is the maximum CO2 uptake rate on the ecosystem scale (𝜇mol CO2 m-2 s-1) (Reichstein et al., 2012). The r is the sum of all the respiratory fluxes in the tower footprint during the sampling period, including respiration from primary producers and microbial communities (Chapin et al., 2006). 

Temperature response curves are a modified Arrhenius equation (Logan, 1982; Falge et al., 2001; Jimenez et al., 2012) (Eq. 2):

NEE ~ a * exp(b*TA)                                     Eq. 2

Where night time NEE is ecosystem respiration (𝜇mol CO2 m-2 s-1), a is the base respiration rate at 0 ℃ (𝜇mol CO2 m-2 s-1), b reveals the sensitivity of NEE to air temperature (TA, ℃). Model parameters are fit using the R package brms (Bürkner, 2017; Padfield et al., 2020). This method allows for prior information on parameter values and uncertainty estimation around predictions and parameters (Bürkner, 2017). Running ~5,000 iterations and four chains with informative priors, we can assess model convergence using posterior predictive checks and Rhat (Bürkner, 2017). 
