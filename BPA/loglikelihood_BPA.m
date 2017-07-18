% Likelihood Approximation Code
% This code calculates the loglikelihood via the BPA method

%% Matrix exponential calculations for the integration
% within household Q matrix calculation
Q=beta_est*Q_bet+gamma_est*Q_gam;
% The within household matrix exponential
eQ=expokit.padm(Q);
%precalculating matrix exponential parts for integration in each household
%and reformatting for efficient integration
eQ1mins(:,:,1)=eQ;
for i=2:numpoints
    eQ1mins(:,:,i)=(eQ)^(1-intmesh(i));
end
eQ1minreshape=reshape(eQ1mins(:,:,end:-1:1),Ssize,numpoints*Ssize);
eQ1permuted=reshape(permute(eQ1mins,[1,3,2]),[Ssize*numpoints,Ssize]);

%distribution of the state of a household after it's first infection
X1_dist=trapz(intmesh,eQ1mins(2,:,:),3);

%% Calculating the Within household likelihood and expected force of infection

%This loop calculates the Expectation for each household while calculating
%each households contribution to the within household likelihood
TotalExphh=zeros(Ht(T_horizon),T_horizon-1);
WithinLHhh=zeros(1,Ht(T_horizon));
for hh=1:Ht(T_horizon)
   [WithinLHhh(hh),TotalExphh(hh,:)]=Epathint(Wmat(hh,tinf(hh):(T_horizon+1)),tinf(hh),Z_mat,Z_mat_cell_col,TotalExphh(hh,:),eQ1permuted,eQ1minreshape,X1_dist,T_horizon,eQ,numpoints,Ssize,bigI,bigintvec);
end
WithinLH=sum(WithinLHhh);
TotalExp=sum(TotalExphh);

%The within household likelihoods for households infected on the last day
for hh=(Ht(T_horizon)+1):Ht(T_horizon+1)
    WithinLH=WithinLH+log(sum(X1_dist.*Z_mat(Wmat(hh,T_horizon+1)+1,:)));
end

%% Calculating Generations of Infected Households

% Distribution of 0th Generation of household infected
h0_dist= [h0col1,poisspdf(cutoffMat,alpha_est*repmat(TotalExp,cutoff,1))];

% The pmf of G calculated via simulation
G_pmf= sim_one_hh([beta_est,gamma_est,alpha_est],simonehhtol);
G_pmf = G_pmf(1:10);

% The convolution matrix
Conv_mat = conv_array(G_pmf,cutoff,cutoff);

% log of the distribution of the number of infected households on each day
ht_dist = log(Conv_mat*h0_dist);

%% Adding together relevant parts for the loglikelihood function

loglikelihood=ht_dist(delHt(2)+1,1);
for t=2:T_horizon
    loglikelihood=loglikelihood+ht_dist(delHt(t+1)+1,t);
end
loglikelihood=loglikelihood+WithinLH;

%%

