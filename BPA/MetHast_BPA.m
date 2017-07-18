%% Metropolis-Hastings Algorithm for BPA Approximation
%The Metropolis-Hastings Algorithm for inferring all parameters of a SIR hh
%model Each simulation is evaluated at different stages
% Inputs:
% HHvec- A vector of number of households to be infected before stopping
% conducting inference.
% currplace- The initial parameter values (or current parameters if adding more
% samples).
% numtogo- The number of MCMC iterations (or the number of iterations left
% if adding more samples)
% Ht_simsave- A vector of the total number of households infected each day from
% simulated data
% Wmat_simsave- A matrix where each columns refer to households and rows
% refer to days. The entries are the number of infections in the household
% on each day.

% Outputs:
% MCMC_samps- The MCMC samples
% runtimes- The run time to get the MCMC samples at each time inference is
% run
function [MCMC_samps,runtimes]=MetHast_BPA(N,currplace,burnin,NUMSAMPLES,Wmat_simsave)
%% Prespecifying the matrix for storing the MCMC samples
MCMC_samps=zeros(3,NUMSAMPLES);

%% Likelihood Calculation Tolerances
% The number of simulations used to calculate the pmf of G
simonehhtol=10^3;
% The number of equally spaced grid points to use for integration
numpoints=25;

%% Prior/ Proposal Distribution

%uniform bounds
lb=[0.05;0.25;0.25];
ub=[1;4;7];

%uniform prior for alpha, R0 and 1/gamma
prior=@(x1,x2,x3)((ub(1)-lb(1))*(lb(1)<x1 & x1<ub(1)))*((ub(2)-lb(2))*(lb(2)<x2 & x2<ub(2)))*((ub(3)-lb(3))*(lb(3)<x3 & x3<ub(3)));
proposalcovmat=diag([0.01,0.02,0.05]);

%Truncated normal proposal distribution
proposal=@(mu,x) mvnpdf(x,mu,proposalcovmat)/mvncdf(lb,ub,mu,proposalcovmat);



%% Extracting some necessary parts from the data for inference
[h,ti]=size(Wmat_simsave);

Ht=zeros(1,ti+1);
for ii=1:ti
    Ht(ii+1)=length(find(Wmat_simsave(:,ii)));
end
Wmat=[zeros(h,1),Wmat_simsave];

[rows,cols]=find(Wmat);
tinf=zeros(1,Ht(end));
for i=1:Ht(end)
    tinf(i)=cols(find(rows==i,1));
end
T_horizon=length(Ht)-1;
delHt=[0,Ht(2:end)-Ht(1:(end-1))];

tic;

%% Initialisation/ Pre-calculations
initialisation_BPA
InitialMCMCval=[currplace(1),currplace(2)/currplace(3),1/currplace(3)]'; %[alpha,R0,1/gamma]
JointMCMCvals=[InitialMCMCval,zeros(3,NUMSAMPLES+burnin)];
alpha_est=currplace(1);
beta_est=currplace(2);
gamma_est=currplace(3);

%% Metropolis-Hastings Algorithm
loglikelihood_BPA
newllh=loglikelihood;
for j=1:(NUMSAMPLES+burnin)
    
    %likelihood values
    oldllh=newllh;
    
    
    %generate candidate within the correct region
    cand_val=[0;0;0];
    while sum(cand_val<=lb) || sum(cand_val>=ub)
        cand_val=mvnrnd(JointMCMCvals(:,j),proposalcovmat)';
    end
    
    %reinterpretting the candidate values in terms of the model
    %parameters
    alpha_est=cand_val(1);
    beta_est=cand_val(2)/cand_val(3);
    gamma_est=1/cand_val(3);
    
    %calculating the candidate likelihood
    loglikelihood_BPA
    candllh=loglikelihood;
    
    %evaluate relevant values of proposal and prior densities
    proposalcgo=proposal(JointMCMCvals(:,j),cand_val);
    proposalogc=proposal(cand_val,JointMCMCvals(:,j));
    priorc=prior(cand_val(1),cand_val(2),cand_val(3));
    prioro=prior(JointMCMCvals(1,j),JointMCMCvals(2,j),JointMCMCvals(3,j));
    
    %acceptance probability
    quantity=exp(candllh-oldllh)*((priorc*proposalogc)/(prioro*proposalcgo));
    u=rand(1,1);
    
    %acceptance rejection step
    if u<quantity
        JointMCMCvals(:,j+1)=cand_val;
        newllh=candllh;
    else
        JointMCMCvals(:,j+1)=JointMCMCvals(:,j);
        newllh=oldllh;
    end
    
end

%saving samples in terms of model parameters
MCMC_samps=[JointMCMCvals(1,(burnin+2):end);JointMCMCvals(2,(burnin+2):end)./JointMCMCvals(3,(burnin+2):end);1./JointMCMCvals(3,(burnin+2):end)];

%saving the run time
runtimes=toc/3600


end
