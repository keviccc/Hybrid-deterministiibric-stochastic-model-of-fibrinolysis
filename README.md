# stochastic
# Hybrid Stochastic Model of Fibrinolysis

A MATLAB implementation of a hybrid deterministic–stochastic model for fibrin clot lysis.

The model describes fibrinolysis at the fibrin-fibre scale by combining deterministic transport/reaction calculations for soluble species with stochastic biochemical reactions on fibrin binding sites. Fibrin degradation is represented through **transverse cutting of fibrin fibres**, allowing the effects of fibre structure and stochastic molecular interactions on clot lysis to be investigated.

## Model overview

Fibrinolysis is governed by interactions between tissue plasminogen activator (tPA), plasminogen (PLG), plasmin (PLS), and fibrin.

Because the number of molecules involved in local reactions on fibrin fibres can be small, these reactions are represented using a stochastic simulation algorithm rather than purely continuous reaction equations.

The model therefore uses a hybrid framework:

- **Deterministic solver** for soluble-phase concentrations
- **Gillespie stochastic simulation** for reactions occurring at fibrin binding sites
- **Discrete fibrin cross-sections** to represent transverse fibre degradation
- Dynamic calculation of:
  - fibrin degradation
  - clot porosity
  - extent of lysis
  - fibrin degradation product distribution
  - bound and free molecular species

A simplified representation of the modelling framework is:

```text
Soluble tPA / PLG / PLS
          |
          | transport and reaction
          v
  Fibrin binding sites
          |
          | stochastic binding,
          | activation and cleavage
          v
 Fibre cross-section degradation
          |
          v
 Porosity + extent of lysis + FDPs
```

## Repository structure

```text
stochastic/
│
├── stochasticTransModel/
│   ├── runTransModel.m
│   ├── funHybridTransModel.m
│   ├── funHybridTransModelExtra.m
│   ├── funHybridTransModelExtra2.m
│   ├── clotPropertiesTrans.m
│   ├── clotPropertiesTransExtra.m
│   └── compareDifferentFibrinolysisPatterns.m
│
├── engineSourceFiles/
│   ├── gillespieSolver.m
│   ├── propensityFunction.m
│   ├── propensityFunctionPar2.m
│   ├── propensityFunctionPar3.m
│   ├── stochiometryMatrix.m
│   ├── stochiometryMatrix2.m
│   ├── hybridConstants.m
│   ├── hybridConstantsExtra.m
│   ├── quasiSteadyPLGCalculation.m
│   ├── rungeKuttaFourth.m
│   ├── rungeKuttaFourth2.m
│   ├── gradient.m
│   └── gradient2.m
│
└── plottingFiles/
    ├── plotAll.m
    └── plotExtra.m
```

### Main files

**`runTransModel.m`**

Main simulation script. Model parameters such as fibrin density, simulation volume, fibre radius, fibre-radius variability, molecular concentrations, simulation time, and lysis threshold can be specified here.

**`funHybridTransModel.m`**

Baseline hybrid fibrinolysis model using fibrin cross-sections of a specified fibre size.

**`funHybridTransModelExtra.m`**

Extended model allowing a distribution of fibrin fibre radii. Fibre sizes are sampled while maintaining the prescribed fibrin density.

**`funHybridTransModelExtra2.m`**

Additional model extension including cellular-volume contributions and associated clot structural properties.

**`gillespieSolver.m`**

Implements the Gillespie stochastic simulation algorithm for fibrin-associated biochemical reactions.

**`propensityFunction.m`**

Defines the stochastic reaction propensities used by the Gillespie solver.

**`stochiometryMatrix2.m`**

Defines changes in molecular populations associated with each stochastic reaction.

**`hybridConstants.m`**

Contains kinetic and structural parameters, including adsorption kinetics, plasminogen activation, fibrinolysis kinetics, fibrin fibre geometry, and initial clot properties.

**`clotPropertiesTrans.m`**

Determines when individual fibrin cross-sections are lysed and calculates quantities including the extent of lysis, voidage, and fibrin degradation product distribution.

## Requirements

The code is written in MATLAB.

Recommended requirements:

- MATLAB
- Statistics and Machine Learning Toolbox for functions such as `normrnd` used for fibre-radius distributions

No compilation is required.

## Running the model

Clone or download the repository and add the model directories to the MATLAB path.

For example:

```matlab
addpath('stochasticTransModel')
addpath('engineSourceFiles')
addpath('plottingFiles')
```

Then run:

```matlab
runTransModel
```

The primary simulation parameters can be changed near the beginning of `runTransModel.m`.

The current example configuration includes:

```matlab
rho_0     = 3;        % fibrin density
U_0       = 100;      % simulation volume
finalTime = 10*60;    % simulation duration [s]
nt        = 1000;     % number of time steps

Cs_tPA = 0.04;
Cs_PLG = 2.2;
Cs_PLS = 0;

lysisLevel = 0.2;

R_f0  = 100;          % mean fibre radius [nm]
sigma = 50;           % fibre-radius standard deviation [nm]
```

The default script currently runs the heterogeneous-fibre model:

```matlab
[t, ploti] = funHybridTransModelExtra( ...
    R_f0, sigma, rho_0, U_0, finalTime, nt, ...
    [Cs_tPA, Cs_PLG, Cs_PLS], ...
    N_fibres, lysisLevel);
```

Because the biochemical reactions are stochastic, repeated simulations with identical parameters may produce slightly different lysis trajectories.

For statistical analyses, multiple realisations should therefore be performed.

## Model outputs

The simulation returns

```matlab
[t, ploti]
```

where `t` contains the simulation times and `ploti` stores the major model outputs.

Important fields include:

| Variable | Description |
|---|---|
| `ploti.CtPA` | Free tPA concentration |
| `ploti.ntPA` | Fibrin-bound tPA |
| `ploti.CPLG` | Free plasminogen concentration |
| `ploti.nPLG` | Fibrin-bound plasminogen |
| `ploti.CPLS` | Free plasmin concentration |
| `ploti.nPLS` | Fibrin-bound plasmin |
| `ploti.LPLS` | Additional fibrin-associated plasmin state |
| `ploti.ntot` | Remaining fibrin binding-site concentration |
| `ploti.epsilon` | Clot voidage / porosity |
| `ploti.extentLysis` | Fractional extent of fibrinolysis |
| `ploti.FDPDist` | Fibrin degradation product distribution |
| `ploti.crossSection` | Evolution of individual fibrin cross-sections |

A basic lysis curve can be plotted using:

```matlab
plot(t, ploti.extentLysis)
xlabel('Time (s)')
ylabel('Extent of lysis')
```

Additional plotting utilities are provided in `plottingFiles/`.

## Hybrid numerical method

During each model time step, two numerical approaches are coupled.

### 1. Continuous-phase calculation

The concentrations of soluble species are updated using a fourth-order Runge–Kutta method.

```text
C(t)
  |
  v
Runge–Kutta solver
  |
  v
C(t + Δt)
```

### 2. Stochastic fibrin reactions

Molecular reactions associated with fibrin cross-sections are simulated using the Gillespie stochastic simulation algorithm.

At each stochastic event:

1. reaction propensities are evaluated;
2. two random numbers are generated;
3. the time to the next reaction is calculated;
4. a reaction is selected according to its relative propensity;
5. molecular populations are updated using the stoichiometric matrix.

This allows discrete molecular events to be retained while the surrounding soluble phase is represented continuously.

## Fibrin degradation

Each fibrin fibre is represented through discrete cross-sections containing fibrin binding sites.

Plasmin-mediated cleavage progressively removes available fibrin sites. Once the number of remaining sites within a cross-section falls below a prescribed threshold, that cross-section is considered lysed.

The loss of cross-sections is subsequently used to calculate:

```text
fibrin degradation
        ↓
extent of lysis
        ↓
change in fibrin volume fraction
        ↓
increase in clot voidage
        ↓
fibrin degradation products
```

This representation is intended to capture the experimentally observed tendency of plasmin-mediated fibrinolysis to produce transverse cleavage of fibrin fibres rather than simple homogeneous radial shrinkage.

## Fibre heterogeneity

`funHybridTransModelExtra.m` allows fibrin fibres to have different radii.

The fibre-radius distribution is controlled by

```matlab
mu
sigma
```

where `mu` is the mean fibre radius and `sigma` determines fibre-size variability.

This allows the effect of heterogeneous fibrin architecture on stochastic lysis behaviour to be investigated while approximately maintaining the prescribed total fibrin density.

## Reproducibility

The stochastic model uses MATLAB's random-number generator.

For reproducible simulations, specify a random seed before running the model:

```matlab
rng(1)
runTransModel
```

For uncertainty or population-level analysis, the model should instead be run repeatedly using independent random-number realisations.

## Scientific background

The model is motivated by experimental and computational studies showing that plasmin-mediated fibrin degradation occurs locally and can result in transverse cutting of fibrin fibres.

Hybrid deterministic–stochastic approaches are particularly useful for fibrinolysis because soluble proteins such as plasminogen may be present at relatively high concentrations, whereas the number of tPA molecules interacting locally with individual fibrin fibres may be small.

Relevant background literature includes:

- Bannish BE, Keener JP, Fogelson AL. *Modelling fibrinolysis: a 3D stochastic multiscale model*. Mathematical Medicine and Biology. 2014;31:17–44.
- Bannish BE et al. *Molecular and Physical Mechanisms of Fibrinolysis and Thrombolysis from Mathematical Modeling and Experiments*. Scientific Reports. 2017;7:6914.

## Notes

This repository contains research code and is primarily intended for mechanistic investigation of fibrinolysis.

Several model variants and analysis scripts are included because the framework has been used to investigate different assumptions regarding fibrin structure, fibre-size distributions, fibrinolysis patterns, and cellular components.

Parameters should therefore be checked carefully before using the code for a new simulation scenario.

## License

No license has currently been specified for this repository.


## Contact

For questions regarding the model or code, please open an issue in this repository.
