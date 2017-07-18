% Some precalculations/ prespecifications that can be run before the
% Metropolis-Hastings algorithm begins

%% Initialisation to allow for Efficient Q-matrix calculations

% Household process sample space size
Ssize=((N+2)*(N+1))/2;

% Matricies Qbet and Qgam are calculated such that Q=beta*Qbet+gamma*Qgam
Q_bet=zeros(Ssize,Ssize);
Q_gam=Q_bet;

% Y is a function that returns the row of the Q matrix corresponding to i
% infectives and r recoverds 
Y=@(i,r) 1+(r.*(2*N+3-r))/2+i;

% I and R give the number of infectives and recovereds in each state
I=zeros(1,Ssize);
R=I;
for i=1:N
    for r=0:(N-i)
        if i~=N
            Q_bet(Y(i,r),Y(i+1,r))=(i*(N-i-r))/(N-1);
        end
        if r~=N
            Q_gam(Y(i,r),Y(i-1,r+1))=i;
        end
        Q_bet(Y(i,r),Y(i,r))=-sum(Q_bet(Y(i,r),:));
        Q_gam(Y(i,r),Y(i,r))=-sum(Q_gam(Y(i,r),:));
        I(Y(i,r))=i;
        R(Y(i,r))=r;
    end
end
R(I==0)=0:N;
p0=[0;1;zeros(Ssize-2,1)];

%% Prespecified matricies of indicators 
%(indicators of which states are possible if n individuals have become infected)
Z_mat=false(N+1,Ssize);
Z_mat_cell_col=cell(1,N+1);
for n=0:N
    Z_mat(n+1,:)=((I+R)==n);
    Z_mat_cell_col{n+1}=repmat(Z_mat(n+1,:)',1,Ssize);
end

%% initialisation for expectation calculations
% grid points for integration in expectation calculations
intmesh=linspace(0,1,numpoints); 

% prespecified matrix for saving matrix exponential parts for integrating
% over
eQ1mins=zeros(Ssize,Ssize,numpoints);

%maximum number of \Delta H_t^{i} considered
cutoff = max(max(delHt)+1,10);

% The distribution of the 0th generation of households on day 1
h0col1=[0;1;zeros(cutoff-2,1)];

% A matrix of the places to evaluate the 0th Generation of infected
% households
cutoffMat=repmat((0:(cutoff-1))',1,T_horizon-1);

% The vector of the number of infected individuals for expectation
% calculation
bigI=repmat(I,1,numpoints);

% Making a vector of Simpsons Integration coffecients
vec_evens_simpson=[];
for i=(2:2:(numpoints-2))
   vec_evens_simpson=[vec_evens_simpson,(1+(i-1)*Ssize):(i*Ssize)]; 
end

vec_odds_simpson=[];
for i=(3:2:(numpoints-1))
   vec_odds_simpson=[vec_odds_simpson,(1+(i-1)*Ssize):(i*Ssize)]; 
end

bigintvec=ones(Ssize*numpoints,1);
bigintvec(vec_evens_simpson)=2;
bigintvec(vec_odds_simpson)=4;


