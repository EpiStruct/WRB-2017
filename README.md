# Walker, Ross and Black, 2017.

This is code in support of the paper: Walker, J., Ross, J.V. and Black A.J. (2017) Inference of epidemiological parameters from household stratified data. https://arxiv.org/abs/1609.09170

If you use or modify this code, please cite the above paper. 
The branching process approximation uses the package Expokit: https://www.maths.uq.edu.au/expokit/
The code expects to find this at ~/MATLAB/+expokit 


## Simulation

A script for simulating household epidemics, SIRhh.m, based on the Ball model simulated by the Doob-Gillespie algorithm.
The simulated data used in the paper is given in the data file, simsforpaper.mat.

## BPA

Perform inference using the branching process approximation (BPA).
The script BPA_example.m loads a test data set, Example_simulation.mat, and runs the Metropolis-Hastings algorithm, MetHast_BPA.m.
 All other scripts and functions in the folder are required to run the Metropolis-Hastings scheme.

## DAMCMC

Perform inference using the Data-augmented MCMC. 
The script DAMCMC_example.m loads a test data set, Example_simulation.mat, and runs the inference routine.
The function DAMCMCfunc.m can be used to run the inference routine on a given FF100 data set.