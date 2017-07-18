# Walker, Ross and Black, 2017.

This is code in support of the paper: Walker, J., Ross, J.V. and Black A.J. (2017) Inference of epidemiological parameters from household stratified data. https://arxiv.org/abs/1609.09170

If you use or modify this code, please cite the above paper. 
The branching process approximation uses the package Expokit: https://www.maths.uq.edu.au/expokit/
The code expects to find this at ~/MATLAB/+expokit 


## simulation

Code for simulation household epidemics. This is based on the Ball model.

## BPA

Perform inference using the branching process approximation (BPA).
The script BPA_example.m loads a test data set and runs the Metropolis-Hastings algorithm.

## DAMCMC

Data-augmented MCMC. 
The script DAMCMC_example.m loads a test data set and runs the inference routine.