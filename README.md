# Continuous Temperature Swing Adsorption Model for CO₂ Capture

This repository contains MATLAB code developed for the simulation and optimization of a continuous Temperature Swing Adsorption (TSA) process for CO₂ capture from biomass combustion flue gas.

The model represents a coupled multi-stage adsorber/desorber system with circulating solid sorbent and steam-assisted regeneration. The main objective is to evaluate how operating parameters such as steam-to-CO₂ ratio, sorbent circulation rate, number of stages, and sorbent properties influence the specific energy demand of the process.

## Thesis Context

**Thesis title:**  
*Continuous Temperature Swing Adsorption for CO₂ Capture from Biomass Combustion Flue Gas: Influence of Source Gas CO₂ Concentration and Sorbent Properties on Process Energy Demand*

The model is based on a steady-state equilibrium-stage approach for continuous counter-current TSA operation.

## Model Overview

The process consists of two connected units:

- Multi-stage adsorber
- Multi-stage desorber
- Circulating solid sorbent loop
- Steam-assisted CO₂ desorption
- Fixed CO₂ capture target
- Energy demand calculation

In the coupled loop, the adsorber outlet loading is passed to the desorber as the rich loading, and the desorber outlet loading is returned to the adsorber as the lean loading.

```text
Adsorber:
lean sorbent in → CO₂ adsorption → rich sorbent out

Desorber:
rich sorbent in → steam regeneration → lean sorbent out
```

## Main Features

- Arbitrary number of adsorber and desorber stages
- Counter-current gas-solid contact
- Toth isotherm implementation
- Coupled adsorber/desorber convergence loop
- Fixed CO₂ capture optimization
- Calculation of required sorbent circulation rate
- Steam-to-CO₂ ratio sweep
- Specific energy demand calculation
- Stagewise CO₂ capture and desorption results
- Stagewise adsorber cooling and desorber heating duties
- Energy demand plots and optimization graphs

## Main Files

| File | Description |
|---|---|
| `run_coupled_TSA_loop.m` | Runs one coupled TSA simulation case |
| `run_energy_optimization.m` | Performs energy optimization by varying steam-to-CO₂ ratio |
| `solve_coupled_TSA_fixed_ms.m` | Solves the coupled TSA loop for fixed sorbent circulation and steam flow |
| `target_sorbent_circulation.m` | Finds the sorbent circulation rate required to reach the target CO₂ capture |
| `solve_adsorber_Nstage.m` | Solves the multi-stage adsorber model |
| `solve_desorber_Nstage.m` | Solves the multi-stage desorber model |
| `Toth_*.m` | Toth isotherm functions and sorbent parameter definitions |

## How to Run

To run a single adsorber or desorber model for any number of stages:

```matlab
run_Nstage_adsorber
run_Nstage_desorber
```

To run a single coupled adsorber/desorber case:

```matlab
run_coupled_TSA_loop
```

To run the energy optimization:

```matlab
run_energy_optimization
```

The optimization script varies the steam-to-CO₂ ratio and calculates the sorbent circulation rate required to maintain the specified CO₂ capture target.

## Optimization Approach

The current optimization fixes the CO₂ capture efficiency and varies the steam-to-CO₂ ratio.

For each steam-to-CO₂ ratio, the required sorbent circulation rate is calculated so that the capture target is achieved.

```text
Input:
steam-to-CO₂ ratio

Constraint:
CO₂ capture efficiency = 90%

Calculated:
required sorbent circulation rate

Output:
specific total heat demand
```

The total specific heat demand includes:

- Desorption heat
- Net sensible heat for sorbent heating/cooling
- Steam generation heat

## Example Result

For the investigated operating range, the minimum specific heat demand was found around:

```text
steam-to-CO₂ ratio         = 0.06 kg steam/kg CO₂
sorbent circulation rate   = 0.076974 kg/s
specific solid circulation = 19.433 kg sorbent/kg CO₂
CO₂ capture efficiency     = 90%
specific heat demand       = 2.962 MJ/kg CO₂
```

This value represents the minimum within the investigated range, not a mathematically proven global minimum.

## Output Plots

The scripts generate plots such as:

- Energy demand vs steam-to-CO₂ ratio
- Energy demand vs specific solid circulation
- 3D optimization curve of steam ratio, solid circulation, and energy demand
- Stagewise CO₂ capture in the adsorber
- Stagewise CO₂ desorption in the desorber
- Stagewise adsorber cooling duty
- Stagewise desorber heating duty

## Notes

The model is a steady-state equilibrium-stage model. It does not currently include detailed hydrodynamic checks, such as minimum fluidization velocity or pressure drop calculations. Therefore, low steam-flow cases may be thermodynamically plausible but still require separate hydrodynamic validation.

## Requirements

The code was developed and tested in MATLAB.

Required MATLAB functionality includes:

- `fsolve`
- `fzero`
- table generation
- basic plotting functions

## Author

Berk Batırbek Çanak

## Academic Use

This repository was developed as part of a master's thesis project on continuous TSA-based CO₂ capture from biomass combustion flue gas.
