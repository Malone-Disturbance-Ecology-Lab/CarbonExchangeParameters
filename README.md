## Carbon Exchange Parameters

Ecosystem CO<sub>2</sub> Flux Partitioning (EC flux partitioning) methods are used to estimate how the total carbon flux between an ecosystem and the atmosphere is divided among different processes, such as photosynthesis (carbon uptake) and respiration (carbon release). These methods are important for understanding the carbon balance of ecosystems. 

Flux partitioning methods can use a light response curve (LRC). We can use daytime NEE (NEE<sub>day</sub>) from tower sites to fit LRCs and nighttime NEE (NEE<sub>night</sub>) to estimate ecosystem respiration with temperature response curves (TRC) at night. 

Here are some common methods for partitioning carbon fluxes:

|Function_Name	|Equation | Data Needed | 
|---------------|-------------------------------------------|---------------|
|LRC_PARMS_01	  | $$NEE = \frac{a_1 \cdot \text{PAR} \cdot a_x}{a_1 \cdot \text{PAR} + a_x} - r $$ | PAR|
|LRC_PARMS_01	  | $$NEE = \frac{aPAR \cdot \text{PAR} \cdot \text{Pmax}}{aPAR \cdot \text{PAR} + \text{Pmax}} - Rd $$ | PAR|
|LRC_PARMS_02	  | $$NEE = \frac{a \cdot \text{PAR} + P_{\text{max}} - \sqrt{(a \cdot \text{PAR} + P_{\text{max}})^2 - 4 \cdot \text{PAR} \cdot a \cdot P_{\text{max}}}}{2 \cdot \Theta} - R_d $$| PAR|
|LRC_PARMS_03	  |$$ NEE = P_\text{max} \left( 1 - \exp \left( -\frac{a \cdot \text{PAR}}{P_{\text{max}} } \right) \right) - R_d $$ | PAR|
|LRC_PARMS_04	  | $$NEE = a \left( 1 - \frac{B \cdot \text{PAR}}{1 + Y \cdot \text{PAR}} \right) - R_d $$| PAR|
|LRC_PARMS_05	  | $$NEE = -a \cdot \exp(-B \cdot \text{PAR}) - Y \cdot \exp(Z \cdot \text{PAR}) $$| PAR|
|LRC_PARMS_06	  |$$NEE = P_m \left( 1 - \exp \left( a \cdot \left( \text{PAR} - I_{\text{comp}} \right) \right) \right) $$ | PAR|
|LRC_PARMS_07	  | $$NEE = \frac{a_1 \cdot \text{PAR} \cdot \beta}{a_1 \cdot \text{PAR} + \beta} $$| PAR; VPD|
