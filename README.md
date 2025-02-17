## Carbon Exchange Parameters

Ecosystem CO<sub>2</sub> Flux Partitioning (EC flux partitioning) methods are used to estimate how the total carbon flux between an ecosystem and the atmosphere is divided among different processes, such as photosynthesis (carbon uptake) and respiration (carbon release). These methods are important for understanding the carbon balance of ecosystems. 

Flux partitioning methods can use a light response curve (LRC). We can use daytime NEE (NEE<sub>day</sub>) from tower sites to fit LRCs and nighttime NEE (NEE<sub>night</sub>) to estimate ecosystem respiration with temperature response curves (TRC) at night. 

Here are some common methods for partitioning carbon fluxes:

|Function_Name	|Equation | Data Needed | 
|---------------|-------------------------------------------|---------------|
|LRC_PARMS_01	  | $$NEE = \frac{a_1 \cdot \text{PAR} \cdot a_x}{a_1 \cdot \text{PAR} + a_x} - r $$ | NEE; PAR|
|LRC_PARMS_01	  | $$NEE = \frac{aPAR \cdot \text{PAR} \cdot \text{Pmax}}{aPAR \cdot \text{PAR} + \text{Pmax}} - Rd $$ | NEE; PAR|
|LRC_PARMS_02	  | $$NEE = \frac{a \cdot \text{PAR} + P_{\text{max}} - \sqrt{(a \cdot \text{PAR} + P_{\text{max}})^2 - 4 \cdot \text{PAR} \cdot a \cdot P_{\text{max}}}}{2 \cdot \Theta} - R_d $$| NEE; PAR|
|LRC_PARMS_03	  | $$NEE = P_{\text{max}} \left( 1 - \exp \left( -\frac{a \cdot \text{PAR}}{P_{\text{max}}} \right) \right) - R_d$$ | NEE; PAR|
|LRC_PARMS_04	  | $$NEE = a \left( 1 - \frac{B \cdot \text{PAR}}{1 + Y \cdot \text{PAR}} \right) - R_d $$| NEE; PAR|
|LRC_PARMS_05	  | $$NEE = -a \cdot \exp(-B \cdot \text{PAR}) - Y \cdot \exp(Z \cdot \text{PAR}) $$| NEE; PAR|
|LRC_PARMS_06	  | $$NEE = P_m \left( 1 - \exp \left( a \cdot \left( \text{PAR} - I_{\text{comp}} \right) \right) \right) $$ | NEE; PAR|
|LRC_PARMS_07	  | $$NEE = \frac{a_1 \cdot \text{PAR} \cdot \beta}{a_1 \cdot \text{PAR} + \beta} $$| NEE; PAR; VPD|
|TRC_PARMS_01   | $$ \text{Reco} = r_b \cdot \exp \left( E_0 \left( \left( \frac{1}{T_{\text{ref}}} - T_0 \right) - \left( \frac{1}{T_{\text{air}}} - T_0 \right) \right) \right)$$| NEE; T|
|TRC_PARMS_02   | $$\text{Reco} = r_0 \cdot \exp \left( \alpha \cdot T + \beta \cdot T^2 \right) $$|NEE; T |
|TRC_PARMS_03   | $$ \text{Reco} = a \left( T - T_{\text{opt}} \right)^2 + E_{\text{Rmax}} $$|NEE; T|
|TRC_PARMS_04   | $$\text{Reco} = R_{\text{ref}} \cdot Q_{10} \cdot \exp \left( \frac{T - T_{\text{ref}}}{10} \right) $$| NEE; T|
|TRC_PARMS_05   | $$\text{Reco} = f_0 + s_1 \cdot \sin(D_{\pi}) + c_1 \cdot \cos(D_{\pi}) + s_2 \cdot \sin(2 \cdot D_{\pi}) + c_2 \cdot \cos(2 \cdot D_{\pi}) $$|NEE; ? |
|TRC_PARMS_06   | $$\text{Reco} = R_{\text{ref}} \cdot \exp \left( \frac{E_a}{R} \left( \frac{1}{T_{\text{ref}}} - \frac{1}{T} \right) \right) $$|NEE; T|

