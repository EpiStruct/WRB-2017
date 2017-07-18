%% Expectation and within household likelihood calculations for a single household
% Inputs:
% W- a vector of the number of infections in the household as of the first
% day of infection
% timehh- the time of the first infection in the household
% Z_mat- Matrix of indicators, row n corresponds to the indicators of possible 
% states if n individuals have been infected
% Z_mat_cell_col- Z_matrix in a format for efficiency of integration% (repeated for each point in the integration for efficiency)
% exp_pathint- Vector for saving the expected path integral of the number of infectious
% individuals of this household on each day
% eQ1permuted-Matrix exponential taken to different powers (for each
% point in the integration) and formated for efficient integration
% eQ1minreshape- Matrix exponential taken to different powers (for each
% point in the integration) and formated for efficient integration
% X1_dist- The distribution of the state of a household after it's first infection
% T_horizon- Time at which inference is being conducted
% eQ- Household matrix exponential
% numpoints- Number of points of integration
% Ssize- Size of the household state space
% bigI- Vector of the number of infected individuals in each state
% (repeated for each point in the integration for efficiency)
% bigintvec- Coefficients for Simpsons integration

% Outputs:
% LHwithin- The within household likelihood associated with this household
% exp_pathint- The expected path integral of the number of infectious
% individuals of this household on each day, conditional upon the within
% household transitions

function [LHwithin,exp_pathint]=Epathint(W,timehh,Z_mat,Z_mat_cell_col,exp_pathint,eQ1permuted,eQ1minreshape,X1_dist,T_horizon,eQ,numpoints,Ssize,bigI,bigintvec)

%% Initial Forward probabilities and Within Household Likelihood Calculation
p_tmin1=X1_dist.*Z_mat(W(1)+1,:);
thenormalisationbit=sum(p_tmin1);
p_tmin1=p_tmin1/thenormalisationbit;
LHwithin=log(thenormalisationbit);

%% Backward Probabilities
sumvec_c=zeros(Ssize,T_horizon+1-timehh);
sumvec_c(:,length(W)-1)=Z_mat(W(end)+1,:)';
for i=(length(W):(-1):3)
    sumvec_c(:,i-2)=Z_mat_cell_col{W(i-1)+1}.*eQ*sumvec_c(:,i-1);
end

%%
for m=1:(T_horizon+1-timehh)
    
    %% Integration
    % Pmf of the state at time s given within household transitions up to
    % time t-1
    feQ=p_tmin1*eQ1minreshape;
    
    % Pmf of the future within household transitions given the state at
    % time s
    eQb=eQ1permuted*sumvec_c(:,m);
    
    % Denominator for Simpsons integration
    denom=feQ((1+(numpoints-1)*Ssize):(numpoints*Ssize))*sumvec_c(:,m);

    % Numerator for Simpsons integration
    numeratorvec=feQ.*bigI.*eQb';
    
    % The integral
    exp_pathint(m+timehh-2)=(numeratorvec*bigintvec)/(denom*3*numpoints);
    
    %% Forward Probability Updates and Within Household Likelihood Calculations
    % Updating forward probability
    p_tmin1=(p_tmin1*eQ).*Z_mat(W(m+1)+1,:);
    thenormalisationbit=sum(p_tmin1);
    p_tmin1=p_tmin1/thenormalisationbit;
    LHwithin=LHwithin+log(thenormalisationbit);
    
end
end