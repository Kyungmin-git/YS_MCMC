# YS_MCMC

MCMC Inversion of Hydrothermal Tremor at Yellowstone Lake

This repository contains MATLAB code for Markov chain Monte Carlo (MCMC) inversion of hydrothermal tremor recorded at Yellowstone Lake.

The inversion uses a gas-pocket forward model to estimate model parameters from the observed tremor spectrum.

## Getting Started

To run the inversion, execute the main MATLAB script:

yellowstone_inversion.m

The repository includes an example tremor spectrum stored in:

YL_Y1701_ELZ_201804_selected_spectra.mat

This dataset can be used as input for testing and reproducing the inversion workflow.

## Code Description

yellowstone_inversion.m: Main script for running the MCMC inversion. 

YL_Y1701_ELZ_201804_selected_spectra.mat: Example tremor spectrum recorded at Yellowstone Lake in April 2018.

forwardmodelgaspocket.m: Computes the theoretical response of the gas-pocket forward model.

compute_G_hankel.m: Computes the radial integration using the Hankel-transform formulation.

generate_excitation.m: Generates MCMC proposals for the excitation parameters.

generatelpparameter.m: Generates MCMC proposals for the physical model parameters.

sigma_propose.m: Generates MCMC proposals for the logarithm of the noise standard deviation.

## Inversion Workflow

The MCMC inversion follows these main steps:

1. Load the observed tremor spectrum.
2. Generate candidate model and excitation parameters.
3. Compute the theoretical spectrum using the gas-pocket forward model and radial integration.
4. Evaluate the agreement between the observed and modeled spectra.
5. Accept or reject candidate parameters according to the MCMC algorithm.
6. Repeat the sampling procedure to characterize the posterior distributions of the model parameters.

## Requirements

- MATLAB
- Any additional MATLAB toolboxes required by the implementation

## Data Availability

An example tremor spectrum is provided in the repository for testing the inversion workflow.


