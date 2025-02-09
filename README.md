## Carbon Exchange Parameters

Flux partitioning methods use a light response curve (LRC). We can use daytime NEE (NEEday) from tower sites to fit LRCs and nighttime NEE (NEEnight) to estimate ecosystem respiration with temperature response curves (TRC). Using daytime data,  fit LRCs:

NEE ~ (a1 * PAR * ax)/(a1 * PAR + ax) - r									
                
Where r is respiration (𝜇mol CO2 m-2 s-1), a1 is the apparent quantum efficiency of CO2 uptake (Reichstein et al., 2012)), PAR is solar radiation (𝜇mol m-2 s-1). ax is the maximum CO2 uptake rate on the ecosystem scale (𝜇mol CO2 m-2 s-1) (Reichstein et al., 2012). The r is the sum of all the respiratory fluxes in the tower footprint during the sampling period, including respiration from primary producers and microbial communities (Chapin et al., 2006). 

Temperature response curves are a modified Arrhenius equation (Logan, 1982; Falge et al., 2001; Jimenez et al., 2012):

NEE ~ a * exp(b*TA)                                   

Where night time NEE is ecosystem respiration (𝜇mol CO2 m-2 s-1), a is the base respiration rate at 0 ℃ (𝜇mol CO2 m-2 s-1), b reveals the sensitivity of NEE to air temperature (TA, ℃). Model parameters are fit using the R package brms (Bürkner, 2017; Padfield et al., 2020). This method allows for prior information on parameter values and uncertainty estimation around predictions and parameters (Bürkner, 2017). Running ~5,000 iterations and four chains with informative priors, we can assess model convergence using posterior predictive checks and Rhat (Bürkner, 2017). 

|Function_Name	|Equation |
|---------------|-------------------------------------------|
|LRC_PARMS_01	  |NEE =  (a1 * PAR * ax)/(a1 * PAR + ax) - r |
|LRC_PARMS_01	  |NEE = (a*PAR*Pmax))/(a*PAR+Pmax)-Rd        |
|LRC_PARMS_02	  |NEE = ((a*PAR+Pmax-(sqrt((a*PAR+Pmax)^2-(4*PAR*a*Pmax))))/(2*Theta))-Rd|
|LRC_PARMS_03	  |NEE = Pmax*(1-exp(-a*PAR/Pmax))-Rd |
|LRC_PARMS_04	  |NEE = a*(1-B*PAR/1+Y*PAR)-Rd|
|LRC_PARMS_05	  |NEE = -a*exp(-B*PAR)-Y*exp(Z*PAR)|
|LRC_PARMS_06	  |NEE = Pm * (1-exp(a*(PAR-Icomp)))|
|LRC_PARMS_07	  |NEE =  (a1 * PAR * beta)/(a1 * PAR + beta) - r  but beta is a function of VPD|
