%An Example of how to run the BPA algorithm for a data set

%load the data set
% Wmat: A matrix with rows corresponding to households and columns
% corresponding to days, entries give the cumulative number of infected individuals have
% infected by the end of the given day in the given household.
load('Example_simulation.mat')

% Number of iterations of burnin and number of samples
burnin=1000;
NUMSAMPLES=1000;

% The household size
N=3;

% Initial parameter guess
currplace=[0.5,0.6,0.4];

% The number of households to become infectious before inference is run

% The Metropolis-Hastings algorithm for the BPA method
[MCMC_samps,runtimes]=MetHast_BPA(N,currplace,burnin,NUMSAMPLES,Wmat);


% Plots of the marginal posterior distributions
burnin=1000;
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
