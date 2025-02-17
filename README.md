## Carbon Exchange Parameters

Ecosystem CO<sub>2</sub> Flux Partitioning (EC flux partitioning) methods are used to estimate how the total carbon flux between an ecosystem and the atmosphere is divided among different processes, such as photosynthesis (carbon uptake) and respiration (carbon release). These methods are important for understanding the carbon balance of ecosystems. 

Flux partitioning methods can use a light response curve (LRC). We can use daytime NEE (NEE<sub>day</sub>) from tower sites to fit LRCs and nighttime NEE (NEE<sub>night</sub>) to estimate ecosystem respiration with temperature response curves (TRC) at night. 

Here are some common methods for partitioning carbon fluxes:

|Function_Name	|Equation |
|---------------|-------------------------------------------|
|LRC_PARMS_01	  | $$NEE = \frac{a_1 \cdot \text{PAR} \cdot a_x}{a_1 \cdot \text{PAR} + a_x} - r $$ |
|LRC_PARMS_01	  |$$NEE = \frac{aPAR \cdot \text{PAR} \cdot Pmax}{aPAR \cdot \text{PAR} + Pmax} - R_d $$        |
|LRC_PARMS_02	  |NEE = ((a*PAR+Pmax-(sqrt((a*PAR+Pmax)^2-(4*PAR*a*Pmax))))/(2*Theta))-Rd|
|LRC_PARMS_03	  |NEE = Pmax*(1-exp(-a*PAR/Pmax))-Rd |
|LRC_PARMS_04	  |NEE = a*(1-B*PAR/1+Y*PAR)-Rd|
|LRC_PARMS_05	  |NEE = -a*exp(-B*PAR)-Y*exp(Z*PAR)|
|LRC_PARMS_06	  |NEE = Pm * (1-exp(a*(PAR-Icomp)))|
|LRC_PARMS_07	  |NEE =  (a1 * PAR * beta)/(a1 * PAR + beta) - r  but beta is a function of VPD|
