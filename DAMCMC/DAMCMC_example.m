% An example of how to run the DA MCMC algorithm

%load the data set 
% Ht: A row vector of the number of newly infected households each day
% Wmat: A matrix with rows corresponding to households and columns
% corresponding to days, entries give the cumulative number of infected individuals have
% infected by the end of the given day in the given household.
load('Example_simulation.mat')
burnin=500;
NUMSAMPLES=1000;

% Household Size
N=3;

% Number of Households
M=50000;

[MCMC_samps,runtime]=DAMCMCfunc(N,M,NUMSAMPLES,burnin,Wmat);

% Plots of the marginal posterior distributions
[al_samps,xal]=ksdensity(MCMC_samps(1,:));
[bet_samps,xbet]=ksdensity(MCMC_samps(2,:));
[gam_samps,xgam]=ksdensity(MCMC_samps(3,:));

subplot(1,3,1)
plot(xal,al_samps)
xlabel('$\alpha$','Interpreter','latex')
subplot(1,3,2)
plot(xbet,bet_samps)
xlabel('$\beta$','Interpreter','latex')
subplot(1,3,3)
plot(xgam,gam_samps)
xlabel('$\gamma$','Interpreter','latex')
