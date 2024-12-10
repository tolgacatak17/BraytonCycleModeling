# Thermodynamic Model of Brayton Cycle (Basic)

This repository contains a project for the **Thermodynamics II (ME466)** course at **Bogazici University**, 
Mechanical Engineering Department. This project contains an engineering model to simulate a Brayton cycle with regeneration and intercooling using various methods in MATLAB.

---

## Problem Description

In this project, a thermodynamic model of a Brayton cycle with enhancements, specifically with regeneration and intercooling, is modeled. The contribution of regeneration and the effects of ambient conditions on specific power output and energetic efficiency are assessed. The schematic of the analyzed system is shown below:

<p align="center">
  <img src="https://github.com/user-attachments/assets/52474981-99f1-466a-b218-2957e48af1e4" alt="System Schematic" width="1000">
</p>

In this thermodynamic model:
- Intercooling and combustion are approximated as heat transfers between the components and the environment. 
- Advanced modeling of these processes, such as liquid water injection and fuel combustion, is completed in the enhanced version of this project.

---

## Assumptions

The following assumptions are made to simplify the modeling:

- Each component of the cycle is analyzed as a control volume at steady state.
- The working fluid is air, and it is treated as an ideal gas.
- Each compressor stage is adiabatic.
- The intercooler exit temperature is equal to T<sub>0</sub> (ambient temperature).
- Intercooling is modeled as heat loss to the surroundings, occurring at the intercooler exit temperature.
- Pressure drop in the regenerator is neglected.
- The regenerator stage is adiabatic.
- Potential and kinetic energy changes are negligible.
- Combustion is modeled as heat transfer to the combustion chamber, occurring at T<sub>max</sub>.
- Heat loss from the turbine occurs at the turbine exit temperature.
- Heat transfer through the regenerator between cold and hot sides occurs at (T<sub>5</sub> + T<sub>8</sub>) / 2.
- Heat loss in the theoretical heat exchanger to the surroundings from state 8 to state 1 occurs at T<sub>1</sub>.

---

## Model Inputs

The base case model inputs are provided below. These values will be varied parametrically for analysis:

<p align="center">
  <img src="https://github.com/user-attachments/assets/4bee0e16-3ab3-42d1-a255-e28702f4962b" alt="Model Inputs" width="700">
</p>


---

## Property Modeling

- **Reference Text**: All property values are based on "*Moran's Principles of Engineering Thermodynamics*".
- **Molecular Weights**: Taken from the table "*Atomic or Molecular Weights and Critical Properties of Selected Elements and Compounds*".
- **Specific Heat (c<sub>p,i</sub>):** Modeled as a function of temperature using the table "*Variation of (c<sub>p</sub>) with Temperature for Selected Ideal Gases*".

- The isentropic turbine efficiency is defined as:

<p align="center">
  <img src="https://github.com/user-attachments/assets/1bdc7144-107c-4cd5-846e-bd6adf5744ab" alt="Isentropic Turbine Efficiency" width="450">
</p>

- **Enthalpy (h) and Entropy (s) Calculations**:

  Enthalpy and entropy are not taken from the tables but are calculated as follows:

  <p align="center">
    <img src="https://github.com/user-attachments/assets/18184d66-b8df-4b70-ab26-5e7d0562e32d" alt="Enthalpy Equation"   width="450">
  </p>

  where h = 0  at T<sub>ref</sub> = 0 K.


  <p align="center">
    <img src="https://github.com/user-attachments/assets/15081e5b-1524-4b42-b6cb-63400917025f" alt="Entropy Equation"   width="450">
  </p>

  where s(T<sub>ref</sub>, P<sub>ref</sub> = 1.70203 kJ/kg.K at T<sub>ref</sub> = 300 K and P<sub>ref</sub> = 1 atm.


- **Constants and Conversions**:

  The constants and conversions are given below:

  <p align="center">
    <img src="https://github.com/user-attachments/assets/22ba788e-e00c-48b5-a209-21bc374f26df" alt="Constants and       Conversions" width="300">
  </p>


---


## Task

The following tasks are to be completed as part of this project:

### Analysis

1. **Base Case**:
   - Analyze the base case conditions and generate state and process tables for verification.

2. **Parametric Analysis 1**:
   - Vary the ambient temperature (T<sub>0</sub>) from 255 K (cold winter night) to 315 K (hot summer afternoon) in 10 K increments, i.e. T<sub>0</sub> = 255, 265, ... , 315 K.
   - Analyze two cases for each temperature:
     1. With regeneration.
     2. Without regeneration (where State 4 and State 5 are identical).
   - Assess the contribution of the regenerator to the cycle performance.

3. **Parametric Analysis 2**:
   - Vary the pressure ratio (r<sub>p,i</sub>) of each compressor stage within a suitable range, (e.g., (r<sub>p,i</sub> = 3, 3.25, 3.5, ...).
   - Analyze two cases for each r<sub>p,i</sub>:
     1. With regeneration.
     2. Without regeneration (where State 4 and State 5 are identical).
   - Assess the contribution of the regenerator to the cycle performance.

### Deliverables

1. **Base Case State and Process Tables**:
   - Produce state and process tables for the base case condition.
   - Verify the accuracy of the model using these tables.

2. **Parametric Analysis 1**:
   - Plot the following metrics against T<sub>0</sub> for cases with and without regeneration:
     - Cycle thermal efficiency (η<sub>thermal</sub>).
     - Backwork ratio (bwr).
     - Specific power output (W<sub>net,specific</sub>).
     - Engine inlet volumetric flow rate (V<sub>in,state1</sub>).
   - Provide **four figures** for this analysis.

3. **Parametric Analysis 2**:
   - Plot the following metrics against r<sub>p,i</sub> for cases with and without regeneration:
     - Cycle thermal efficiency (η<sub>thermal</sub>).
     - Backwork ratio (bwr).
     - Specific power output (W<sub>net,specific</sub>).
     - Engine inlet volumetric flow rate (V<sub>in,state1</sub>).
   - Provide **four figures** for this analysis.

### Discussion

- What is the typical backwork ratio (bwr) for gas turbines? Discuss the values obtained in this project.
- What is the suitable range of \(r<sub>p,i</sub>\) for compressor stages? Explain how this range was determined.
- Discuss the effects of the regenerator on bwr, specific power output, and thermal efficiency.
