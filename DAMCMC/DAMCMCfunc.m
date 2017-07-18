%% Data Augmented Method for SIR Household Model
% Inputs:
% sampo- The sample index
% N- The household size
% M- The number of households
% hh_ind- The index of the number of infected households we are considering
% NUMSAMPLES- The number of MCMC samples collected
% burnin- The number of iterations of burn-in

% Outputs:
% mcmc_samps- The MCMC samples
% runtime- The run time of the DAMCMC algorithm
function [mcmc_samps,runtime]=DAMCMCfunc(N,M,NUMSAMPLES,burnin,Wmat)

%% Some Initialisation/ Prespecification

%Proposal type pmf
movechoose=[0.05,0.05,0.3,0.3,0.3];

%Vector to keep track of each kind of acceptance in the Hastings step
no_of_acceptances=zeros(1,5);

% Extracting the relevant parts of the data files
[h,ti]=size(Wmat);
Ht=zeros(1,ti);
for ii=1:ti
    Ht(ii)=length(find(Wmat(:,ii)));
end

delWmat=[Wmat(:,1)';(Wmat(:,2:end)-Wmat(:,1:(end-1)))'];
T_horizon=length(Ht);
hhnumbo=Ht(T_horizon);
data=delWmat(1:T_horizon,1:hhnumbo);

% the number of iterations between display updates
update_its=200000;
mcmc_samps=zeros(NUMSAMPLES,3);

%% Some General Data Manipulation/Initialisation
totalI=sum(sum(data));
[eventday,eventhh]=find(data==1);
[eventday2,eventhh2]=find(data==2);
[eventday3,eventhh3]=find(data==3);
data_mod_temp=[eventday,eventhh;eventday2,eventhh2;eventday2,eventhh2;eventday3,eventhh3;eventday3,eventhh3;eventday3,eventhh3];
[~,sorto]=sort(data_mod_temp(:,1));
data_mod=[data_mod_temp(sorto,:),ones(totalI,1)];
for i=1:hhnumbo
    data_mod(find(data_mod(:,2)==i,1),3)=2;
end
%% Prior/ Proposal Distributions for Gamma and Beta

%Beta profile likelihood
lhprop_bet=@(pathint,tot_w_I) gamrnd(tot_w_I+1,1/pathint);

%Gamma profile likelihood
lhprop_gam=@(pathint,numrecs) gamrnd(numrecs+1,1/pathint);

%Alph profile likelihood
lhprop_alph=@(pathint,tot_b_I) gamrnd(tot_b_I,1/pathint);

%Proposal Distribution for First Infection time
I1prop=@(bet,alph,gam,I2) log(rand(1,1)*(exp((bet+gam+alph)*min(1,I2))-1)+1)/(bet+gam+alph);

%The Augmented Log Likelihood
logIR_jointdist=@(bet,alph,gam,path_w_inf,path_b_inf,pathrec,wmix,bmix,recmix) sum(log((bet/(N-1))*wmix))+sum(log(alph*bmix./(N*(M-1))))+sum(log(gam*recmix))-path_w_inf-path_b_inf-pathrec;


%Initial Beta and Gamma values
bet_old=0.3;
alph_old=0.2 ;
gam_old=0.2;


%% Initialisation of missing data

%Makes sure the epidemic doesn't die out
numinf=0;
while sum(numinf<=0) 
    
    %Event times vector
    IR=zeros(1,2*totalI);
    
    %Tracking the households in which events occur
    the_hh_with_transition=zeros(1,2*totalI);
    
    %Tracking the transition types: 2=between household transition, 1=within
    %household transition, 3=recovery transition
    transtype=zeros(1,2*totalI);
    
    %Index corresponding to the number of transitions already considered+1
    countthetrans=1;
    
    % This generates all infection and recovery transitions in each
    % household (we initialise such that only the first transition is
    % bet_hh and all recoveries occur at some exp distributed time)
    for hh=1:hhnumbo
        
        % Loop to ensure feasibility within each household
        impossible=1;
        while impossible==1
            %Vector of the days on which infections occur
            timesforhh=data_mod(data_mod(:,2)==hh,1);
            
            % The infection times
            IR(countthetrans:(countthetrans+length(timesforhh)-1))=sort(timesforhh'-1+rand(1,length(timesforhh)));
            
            % The recovery times
            IR((countthetrans+length(timesforhh)):(countthetrans+2*length(timesforhh)-1))=IR(countthetrans:(countthetrans+length(timesforhh)-1))+exprnd(1/gam_old,1,length(timesforhh));
            
            % The transition types (we begin by assuming the first
            % infection is between household and the other infections are
            % within household)
            transtype(countthetrans:(countthetrans+2*length(timesforhh)-1))=[2,ones(1,length(timesforhh)-1),3*ones(1,length(timesforhh))];
            
            % Tracking the household in which these transitions occur
            the_hh_with_transition(countthetrans:(countthetrans+2*length(timesforhh)-1))=hh;
            
            % Checking feasibility within the household
            [~,shuffler]=sort(IR(countthetrans:(countthetrans+2*length(timesforhh)-1)));
            trans_temp=transtype(countthetrans:(countthetrans+2*length(timesforhh)-1));
            trans_temp=trans_temp(shuffler);
            inf_in_hh=cumsum((trans_temp~=3)-(trans_temp==3));
            if ~sum(inf_in_hh(1:end-1)==0)
                impossible=0;
            end
        end
        % Updating the index
        countthetrans=countthetrans+2*length(timesforhh);
    end
    
    %Reordering the transitions
    [IR,shuffler]=sort(IR);   
    transtype=transtype(shuffler);
    the_hh_with_transition=the_hh_with_transition(shuffler);
    
    % Only save tranisitions that occur before the epidemic runs out
    considero=sum(IR<T_horizon);
    IR=IR(1:considero);
    transtype=transtype(1:considero);
    the_hh_with_transition=the_hh_with_transition(1:considero);
    
    %Total number of infectives at each transition time
    numinf=cumsum((transtype~=3)-(transtype==3))'; 
end
%interarrival times until T_horizon
timeintervals=([IR(2:end),T_horizon]-IR);

%Number of susceptibles and infectives in each of the infected households at each transition 
sushh=N*ones(length(numinf),hhnumbo);
infhh=zeros(length(numinf),hhnumbo);
for ind_trans=1:length(numinf)
    infhh(ind_trans:end,the_hh_with_transition(ind_trans))=infhh(ind_trans,the_hh_with_transition(ind_trans))+(transtype(ind_trans)==1)+(transtype(ind_trans)==2)-(transtype(ind_trans)==3);
    sushh(ind_trans:end,the_hh_with_transition(ind_trans))=sushh(ind_trans,the_hh_with_transition(ind_trans))-(transtype(ind_trans)==1)-(transtype(ind_trans)==2);
end
I_matto=repmat(numinf,1,hhnumbo)-infhh;

% The path integrals of the rate for each transition type (for the exponential term in
% the Augmented likelihood)
path_w_infbit=bet_old*timeintervals*(sum(infhh.*sushh,2)/(N-1));
path_b_infbit=alph_old*timeintervals*((numinf*N*(M-hhnumbo)+sum(I_matto.*sushh,2))/(N*(M-1)));
pathrecbit=gam_old*timeintervals*numinf;

% Terms for the cofficients in the product of the likelihood expression: si, s(I-i) and i.
modtrans=transtype(2:end);
hhmod=the_hh_with_transition(2:end);
wmix=sushh(find(modtrans==1)+(hhmod(modtrans==1)-1)*length(transtype)).*infhh(find(modtrans==1)+(hhmod(modtrans==1)-1)*length(transtype));
bmix=sushh(find(modtrans==2)+(hhmod(modtrans==2)-1)*length(transtype)).*(numinf(modtrans==2)'-infhh(find(modtrans==2)+(hhmod(modtrans==2)-1)*length(transtype)));
rmix=infhh(find(modtrans==3)+(hhmod(modtrans==3)-1)*length(transtype));

%% Gibbs Sampling/ Single Component MH
tic;
for jj=2:(NUMSAMPLES+burnin)
    
    
    %% GIBBS SAMPLES
    
    %Gibbs Sampling for Alpha
    alph_new=lhprop_alph(path_b_infbit/alph_old,sum(transtype==2));
    %while alph_new>2 || alph_new<0
    while alph_new>1 || alph_new<0.05
        alph_new=lhprop_alph(path_b_infbit/alph_old,sum(transtype==2));
    end
    path_b_infbit=(alph_new/alph_old)*path_b_infbit;
    alph_old=alph_new;
    
    %Gibbs for Gamma
    gam_new=lhprop_gam(pathrecbit/gam_old,sum(transtype==3)-2);
    while 1/gam_new>7 || 1/gam_new<0.25
        gam_new=lhprop_gam(pathrecbit/gam_old,sum(transtype==3)-2);
    end
    pathrecbit=(gam_new/gam_old)*pathrecbit;
    gam_old=gam_new;
    
    
    %Gibbs Sampling for Beta
    bet_new=lhprop_bet(path_w_infbit/bet_old,sum(transtype==1));
    while bet_new/gam_new>4 || bet_new/gam_new<0.25
        bet_new=lhprop_bet(path_w_infbit/bet_old,sum(transtype==1));
    end
    path_w_infbit=(bet_new/bet_old)*path_w_infbit;
    bet_old=bet_new;
    
    %Gibbs for I1 and some recalculating
    IR(1)=I1prop(bet_new,alph_new,gam_new,IR(2));
    temp_int1=timeintervals(1);
    timeintervals(1)=IR(2)-IR(1);
    path_w_infbit=path_w_infbit+(timeintervals(1)-temp_int1)*bet_new;
    path_b_infbit=path_b_infbit+(timeintervals(1)-temp_int1)*alph_new;
    pathrecbit=pathrecbit+(timeintervals(1)-temp_int1)*gam_new;
    
    logold_IR_dens=logIR_jointdist(bet_new,alph_new,gam_new,path_w_infbit,path_b_infbit,pathrecbit,wmix,bmix,rmix);
    
    
    %% Single Component MH
    hastingsstep=rand(1,1);
    if hastingsstep<movechoose(1)
        %% Infection Shift
        IR_temp=IR;
        
        %randomly choose an infection event (not the first infection)
        whichI=1+ceil(rand(1,1)*(totalI-1)); 
        indexofI=find(cumsum(transtype~=3)==whichI,1);
        
        %uniformly shift the infection over the day
        if data_mod(whichI)~=1
            IR_temp(indexofI)=data_mod(whichI)-1+rand(1,1); %choose U(0,min(1,rectime)) candidate on relevant day
        else
            IR_temp(indexofI)=IR_temp(1)+rand(1,1)*min(1-IR_temp(1)); %choose U(I1,min(1,rectime)) candidate on relevant day
        end
        
        %reorder the transition times etc. as necessary
        [IR_temp,shuffler]=sort(IR_temp);
        transtype_temp=transtype(shuffler);
        the_hh_with_transition_temp=the_hh_with_transition(shuffler);
        
        %%%check feasibility in the household
        impossible=0;
        %find all transitions for this household
        hh_transo=transtype_temp(the_hh_with_transition_temp==the_hh_with_transition(indexofI));
        %the first entry must be a between household transition
        if hh_transo(1)~=2 
            impossible=1;
        else
            %calculating the total number of infectious individuals at each
            %transition event
            inf_in_hh=cumsum((hh_transo~=3)-(hh_transo==3));
            %checks that if a household dies out that the next transition
            %is a between household transition
            if sum(hh_transo(logical([0,inf_in_hh(1:end-1)==0]))~=2)
                impossible=1;
            end
        end
        
        %check the epidemic doesn't die out
        numinf_temp=cumsum((transtype_temp~=3)-(transtype_temp==3))';
        if ~sum(numinf_temp<=0) && impossible==0
            
            %recalculate the states of each household at each transition
            indextranso=find(shuffler==indexofI,1);
            if indextranso>indexofI
                infhh_temp=[infhh(1:(indexofI-1),:);infhh((indexofI+1):indextranso,:);infhh(indextranso:end,:)];
                infhh_temp(indexofI:(indextranso-1),the_hh_with_transition(indexofI))=infhh_temp(indexofI:(indextranso-1),the_hh_with_transition(indexofI))-1;
                
                sushh_temp=[sushh(1:(indexofI-1),:);sushh((indexofI+1):indextranso,:);sushh(indextranso:end,:)];
                sushh_temp(indexofI:(indextranso-1),the_hh_with_transition(indexofI))=sushh_temp(indexofI:(indextranso-1),the_hh_with_transition(indexofI))+1;
                
            elseif indextranso<indexofI
                infhh_temp=[infhh(1:(indextranso-1),:);infhh((indextranso-1):(indexofI-1),:);infhh((indexofI+1):end,:)];
                infhh_temp(indextranso:indexofI,the_hh_with_transition(indexofI))=infhh_temp(indextranso:indexofI,the_hh_with_transition(indexofI))+1;
                
                sushh_temp=[sushh(1:(indextranso-1),:);sushh((indextranso-1):(indexofI-1),:);sushh((indexofI+1):end,:)];
                sushh_temp(indextranso:indexofI,the_hh_with_transition(indexofI))=sushh_temp(indextranso:indexofI,the_hh_with_transition(indexofI))-1;
            else
                infhh_temp=infhh;
                sushh_temp=sushh;
            end
            I_matto_temp=repmat(numinf_temp,1,hhnumbo)-infhh_temp;
            
            
            
            %recalculate path integrals and interarrival times
            timeintervals_temp=([IR_temp(2:end),T_horizon]-IR_temp);
            path_w_infbit_temp=bet_new*timeintervals_temp*(sum(infhh_temp.*sushh_temp,2)/(N-1));
            path_b_infbit_temp=alph_new*timeintervals_temp*((numinf_temp*N*(M-hhnumbo)+sum(I_matto_temp.*sushh_temp,2))/(N*(M-1)));
            pathrecbit_temp=gam_new*timeintervals_temp*numinf_temp;
            
            %recalculate likelihood coefficient terms
            modtrans=transtype_temp(2:end);
            hhmod=the_hh_with_transition_temp(2:end);
            wmix_temp=sushh_temp(find(modtrans==1)+(hhmod(modtrans==1)-1)*length(transtype_temp)).*infhh_temp(find(modtrans==1)+(hhmod(modtrans==1)-1)*length(transtype_temp));
            bmix_temp=sushh_temp(find(modtrans==2)+(hhmod(modtrans==2)-1)*length(transtype_temp)).*(numinf_temp(modtrans==2)'-infhh_temp(find(modtrans==2)+(hhmod(modtrans==2)-1)*length(transtype_temp)));
            rmix_temp=infhh_temp(find(modtrans==3)+(hhmod(modtrans==3)-1)*length(transtype_temp));
            
            %calculate the candidate augmented likelihood
            acceptance=logIR_jointdist(bet_new,alph_new,gam_new,path_w_infbit_temp,path_b_infbit_temp,pathrecbit_temp,wmix_temp,bmix_temp,rmix_temp);
            
            %acception rejection steps
            r=rand(1,1);
            if r<exp(acceptance-logold_IR_dens)
                
                %%%saving candidate values
                
                %tracking the acceptances
                no_of_acceptances(1)=no_of_acceptances(1)+1;
                
                %augmented data
                IR=IR_temp;
                transtype=transtype_temp;
                the_hh_with_transition=the_hh_with_transition_temp;
                
                %path integrals
                path_b_infbit=path_b_infbit_temp;
                path_w_infbit=path_w_infbit_temp;
                pathrecbit=pathrecbit_temp;
                timeintervals=timeintervals_temp;
                
                %likelihood coefficient terms
                rmix=rmix_temp;
                wmix=wmix_temp;
                bmix=bmix_temp;
                
                %states at each transition
                numinf=numinf_temp;
                sushh=sushh_temp;
                infhh=infhh_temp;
            end
        end
        
    elseif hastingsstep<sum(movechoose(1:2))
        
        %% Infection type swap
        transtype_temp=transtype;
        
        %randomly choose an infection event
        whichI=1+ceil(rand(1,1)*(totalI-1));
        indexofI=find(cumsum(transtype~=3)==whichI,1);
        
        %swap its type
        if transtype(indexofI)==2
            transtype_temp(indexofI)=1;
        elseif transtype(indexofI)==1
            transtype_temp(indexofI)=2;
        end
        
        %check feasibility of household transitions
        impossible=0;
        hh_transo=transtype_temp(the_hh_with_transition==the_hh_with_transition(indexofI));
        if hh_transo(1)~=2
            impossible=1;
        else
            inf_in_hh=cumsum((hh_transo~=3)-(hh_transo==3));
            if sum(hh_transo(logical([0,inf_in_hh(1:end-1)==0]))~=2)
                impossible=1;
            end
        end
        
        
        if impossible==0
            
            %recalculate likelihood coefficient terms
            modtrans=transtype_temp(2:end);
            hhmod=the_hh_with_transition(2:end);
            wmix_temp=sushh(find(modtrans==1)+(hhmod(modtrans==1)-1)*length(transtype_temp)).*infhh(find(modtrans==1)+(hhmod(modtrans==1)-1)*length(transtype_temp));
            bmix_temp=sushh(find(modtrans==2)+(hhmod(modtrans==2)-1)*length(transtype_temp)).*(numinf(modtrans==2)'-infhh(find(modtrans==2)+(hhmod(modtrans==2)-1)*length(transtype_temp)));
            
            %calculate the candidate loglikelihood
            acceptance=logIR_jointdist(bet_new,alph_new,gam_new,path_w_infbit,path_b_infbit,pathrecbit,wmix_temp,bmix_temp,rmix);
            
            %acceptance/rejection step
            r=rand(1,1);
            if r<exp(acceptance-logold_IR_dens)
                
                %%% saving candidate values
                
                %tracking the acceptances
                no_of_acceptances(2)=no_of_acceptances(2)+1;
                
                %transitions
                transtype=transtype_temp;
                
                %likelihood coefficient terms
                wmix=wmix_temp;
                bmix=bmix_temp;
                
            end
        end
        
    elseif hastingsstep<sum(movechoose(1:3))
        %% Recovery Shift
        IR_temp=IR;
        
        %choose a recovery event
        whichR=ceil(rand(1,1)*(sum(transtype==3)));
        indexofR=find(cumsum(transtype==3)==whichR,1);
        
        %find the time of the first event in the household
        hh_time=IR(find(the_hh_with_transition==the_hh_with_transition(indexofR),1));
        
        %generate a uniformly distributed recovery candidate
        IR_temp(indexofR)=hh_time+rand(1,1)*(T_horizon-hh_time);
        
        %reorder transition times as necessary
        [IR_temp,shuffler]=sort(IR_temp);
        transtype_temp=transtype(shuffler);
        the_hh_with_transition_temp=the_hh_with_transition(shuffler);
        
        %check household feasibility
        impossible=0;
        hh_transo=transtype_temp(the_hh_with_transition_temp==the_hh_with_transition(indexofR));
        if hh_transo(1)~=2
            impossible=1;
        else
            inf_in_hh=cumsum((hh_transo~=3)-(hh_transo==3));
            if sum(hh_transo(logical([0,inf_in_hh(1:end-1)==0]))~=2) 
                impossible=1;
            end
        end
        
        %check the epidemic doesn't die out
        numinf_temp=cumsum((transtype_temp~=3)-(transtype_temp==3))';
        if ~sum(numinf_temp<=0) && impossible==0
            
            %recalculate the states of each household at each transition
            indextranso=find(shuffler==indexofR,1);
            if indextranso>indexofR
                infhh_temp=[infhh(1:(indexofR-1),:);infhh((indexofR+1):indextranso,:);infhh(indextranso:end,:)];
                infhh_temp(indexofR:(indextranso-1),the_hh_with_transition(indexofR))=infhh_temp(indexofR:(indextranso-1),the_hh_with_transition(indexofR))+1;
                
                sushh_temp=[sushh(1:(indexofR-1),:);sushh((indexofR+1):indextranso,:);sushh(indextranso:end,:)];
                
            elseif indextranso<indexofR
                infhh_temp=[infhh(1:(indextranso-1),:);infhh((indextranso-1):(indexofR-1),:);infhh((indexofR+1):end,:)];
                infhh_temp(indextranso:indexofR,the_hh_with_transition(indexofR))=infhh_temp(indextranso:indexofR,the_hh_with_transition(indexofR))-1;
                
                sushh_temp=[sushh(1:(indextranso-1),:);sushh((indextranso-1):(indexofR-1),:);sushh((indexofR+1):end,:)];
            else
                infhh_temp=infhh;
                sushh_temp=sushh;
            end
            I_matto_temp=repmat(numinf_temp,1,hhnumbo)-infhh_temp;
            
            %recalculate path integrals
            timeintervals_temp=([IR_temp(2:end),T_horizon]-IR_temp);
            path_w_infbit_temp=bet_new*timeintervals_temp*(sum(infhh_temp.*sushh_temp,2)/(N-1));
            path_b_infbit_temp=alph_new*timeintervals_temp*((numinf_temp*N*(M-hhnumbo)+sum(I_matto_temp.*sushh_temp,2))/(N*(M-1)));
            pathrecbit_temp=gam_new*timeintervals_temp*numinf_temp;
            
            %recalculate likelihood coefficient terms
            modtrans=transtype_temp(2:end);
            hhmod=the_hh_with_transition_temp(2:end);
            wmix_temp=sushh_temp(find(modtrans==1)+(hhmod(modtrans==1)-1)*length(transtype_temp)).*infhh_temp(find(modtrans==1)+(hhmod(modtrans==1)-1)*length(transtype_temp));
            bmix_temp=sushh_temp(find(modtrans==2)+(hhmod(modtrans==2)-1)*length(transtype_temp)).*(numinf_temp(modtrans==2)'-infhh_temp(find(modtrans==2)+(hhmod(modtrans==2)-1)*length(transtype_temp)));
            rmix_temp=infhh_temp(find(modtrans==3)+(hhmod(modtrans==3)-1)*length(transtype_temp));
            
            %calculuate the loglikelihood
            acceptance=logIR_jointdist(bet_new,alph_new,gam_new,path_w_infbit_temp,path_b_infbit_temp,pathrecbit_temp,wmix_temp,bmix_temp,rmix_temp);
            
            %acceptance rejection step
            r=rand(1,1);
            if r<exp(acceptance-logold_IR_dens)
                
                %%% saving candidate values
                
                %tracking the acceptances
                no_of_acceptances(3)=no_of_acceptances(3)+1;
                
                %augmented data
                IR=IR_temp;
                transtype=transtype_temp;
                the_hh_with_transition=the_hh_with_transition_temp;
                
                %path integrals
                path_b_infbit=path_b_infbit_temp;
                path_w_infbit=path_w_infbit_temp;
                pathrecbit=pathrecbit_temp;
                timeintervals=timeintervals_temp;
                
                %likelihood coefficient terms
                rmix=rmix_temp;
                wmix=wmix_temp;
                bmix=bmix_temp;
                
                %states at each transition
                numinf=numinf_temp;
                sushh=sushh_temp;
                infhh=infhh_temp;
            end
        end
        
        
    elseif hastingsstep<sum(movechoose(1:4))
        
        %% Recovery Removal
        if sum(transtype==3)~=0
            IR_temp=IR;
            
            %randomly select a recovery event
            whichR=ceil(rand(1,1)*(sum(transtype==3)));
            indexofR=find(cumsum(transtype==3)==whichR,1);
            hh_time=IR(find(the_hh_with_transition==the_hh_with_transition(indexofR),1));
            
            %remove the recovery
            IR_temp(indexofR)=[];
            transtype_temp=transtype;
            transtype_temp(indexofR)=[];
            the_hh_with_transition_temp=the_hh_with_transition;
            the_hh_with_transition_temp(indexofR)=[];
            
            %recalculate the states of each household at each transition
            numinf_temp=cumsum((transtype_temp~=3)-(transtype_temp==3))';
            sushh_temp=[sushh(1:(indexofR-1),:);sushh(indexofR+1:end,:)];
            infhh_temp=[infhh(1:(indexofR-1),:);infhh(indexofR+1:end,:)];
            infhh_temp(indexofR:end,the_hh_with_transition(indexofR))=infhh_temp(indexofR:end,the_hh_with_transition(indexofR))+1;
            I_matto_temp=repmat(numinf_temp,1,hhnumbo)-infhh_temp;
            
            %recalculate path integrals
            timeintervals_temp=([IR_temp(2:end),T_horizon]-IR_temp);
            path_w_infbit_temp=bet_new*timeintervals_temp*(sum(infhh_temp.*sushh_temp,2)/(N-1)); %path integral related to infections
            path_b_infbit_temp=alph_new*timeintervals_temp*((numinf_temp*N*(M-hhnumbo)+sum(I_matto_temp.*sushh_temp,2))/(N*(M-1)));
            pathrecbit_temp=gam_new*timeintervals_temp*numinf_temp; %path integral related to recoveries
            
            %recalculate the likelihood coefficient terms
            modtrans=transtype_temp(2:end);
            hhmod=the_hh_with_transition_temp(2:end);
            wmix_temp=sushh_temp(find(modtrans==1)+(hhmod(modtrans==1)-1)*length(transtype_temp)).*infhh_temp(find(modtrans==1)+(hhmod(modtrans==1)-1)*length(transtype_temp));
            bmix_temp=sushh_temp(find(modtrans==2)+(hhmod(modtrans==2)-1)*length(transtype_temp)).*(numinf_temp(modtrans==2)'-infhh_temp(find(modtrans==2)+(hhmod(modtrans==2)-1)*length(transtype_temp)));
            rmix_temp=infhh_temp(find(modtrans==3)+(hhmod(modtrans==3)-1)*length(transtype_temp));
            
            %calculate the candidate likelihood
            acceptance=logIR_jointdist(bet_new,alph_new,gam_new,path_w_infbit_temp,path_b_infbit_temp,pathrecbit_temp,wmix_temp,bmix_temp,rmix_temp);
            
            %calculate the ratio of proposal densities
            propratio=(sum(transtype==3)*movechoose(5))/((T_horizon-hh_time)*hhnumbo*movechoose(4));
            
            %acceptance/rejection step
            r=rand(1,1);
            if r<exp(acceptance-logold_IR_dens)*propratio
                
                %%% save candidate values
                
                %tracking the acceptances
                no_of_acceptances(4)=no_of_acceptances(4)+1;
                
                %augmented data
                IR=IR_temp;
                transtype=transtype_temp;
                the_hh_with_transition=the_hh_with_transition_temp;
                
                %path integrals
                path_b_infbit=path_b_infbit_temp;
                path_w_infbit=path_w_infbit_temp;
                pathrecbit=pathrecbit_temp;
                timeintervals=timeintervals_temp;
                
                %likelihood coefficients
                rmix=rmix_temp;
                wmix=wmix_temp;
                bmix=bmix_temp;
                
                %states at each transition
                numinf=numinf_temp;
                sushh=sushh_temp;
                infhh=infhh_temp;
            end
        end
        
    else
        
        %% Recovery Insertion
        if sum(transtype==3)<totalI
            
            %randomly select a household
            the_intro_hh=ceil(rand(1,1)*hhnumbo);
            
            %insert a uniformly distributed recovery time
            hh_time=IR(find(the_hh_with_transition==the_intro_hh,1));
            IR_temp=[IR,hh_time+rand(1,1)*(T_horizon-hh_time)];
            the_hh_with_transition_temp=[the_hh_with_transition,the_intro_hh];
            transtype_temp=[transtype,3];
            
            %reorder transitions
            [IR_temp,shuffler]=sort(IR_temp);
            transtype_temp=transtype_temp(shuffler);
            the_hh_with_transition_temp=the_hh_with_transition_temp(shuffler);
            
            
            %check household feasibility
            impossible=0;
            hh_transo=transtype_temp(the_hh_with_transition_temp==the_intro_hh);
            if hh_transo(1)~=2
                impossible=1;
            else
                inf_in_hh=cumsum((hh_transo~=3)-(hh_transo==3));
                if sum(hh_transo(logical([0,inf_in_hh(1:end-1)==0]))~=2) 
                    impossible=1;
                end
            end
            
            %check the epidemic doesn't die out
            numinf_temp=cumsum((transtype_temp~=3)-(transtype_temp==3))';
            if ~sum(numinf_temp<=0) && impossible==0
                
                %recalculate the states of each household at each transition
                indextranso=find(shuffler==length(numinf_temp),1);
                infhh_temp=[infhh(1:(indextranso-1),:);infhh(indextranso-1,:);infhh(indextranso:end,:)];
                sushh_temp=[sushh(1:(indextranso-1),:);sushh(indextranso-1,:);sushh(indextranso:end,:)];
                infhh_temp(indextranso:end,the_intro_hh)=infhh_temp(indextranso:end,the_intro_hh)-1;
                
                I_matto_temp=repmat(numinf_temp,1,hhnumbo)-infhh_temp;
                
                %recalculate path integrals
                timeintervals_temp=([IR_temp(2:end),T_horizon]-IR_temp);
                path_w_infbit_temp=bet_new*timeintervals_temp*(sum(infhh_temp.*sushh_temp,2)/(N-1)); %path integral related to infections
                path_b_infbit_temp=alph_new*timeintervals_temp*((numinf_temp*N*(M-hhnumbo)+sum(I_matto_temp.*sushh_temp,2))/(N*(M-1)));
                pathrecbit_temp=gam_new*timeintervals_temp*numinf_temp; %path integral related to recoveries
                
                %recalculate likelihood coefficient terms
                modtrans=transtype_temp(2:end);
                hhmod=the_hh_with_transition_temp(2:end);
                wmix_temp=sushh_temp(find(modtrans==1)+(hhmod(modtrans==1)-1)*length(transtype_temp)).*infhh_temp(find(modtrans==1)+(hhmod(modtrans==1)-1)*length(transtype_temp));
                bmix_temp=sushh_temp(find(modtrans==2)+(hhmod(modtrans==2)-1)*length(transtype_temp)).*(numinf_temp(modtrans==2)'-infhh_temp(find(modtrans==2)+(hhmod(modtrans==2)-1)*length(transtype_temp)));
                rmix_temp=infhh_temp(find(modtrans==3)+(hhmod(modtrans==3)-1)*length(transtype_temp));
                
                %calculate the candidate loglikelihood
                acceptance=logIR_jointdist(bet_new,alph_new,gam_new,path_w_infbit_temp,path_b_infbit_temp,pathrecbit_temp,wmix_temp,bmix_temp,rmix_temp);
                
                %calculate the ratio of proposal densities
                propratio=((T_horizon-hh_time)*hhnumbo*movechoose(4))/((sum(transtype==3)+1)*movechoose(5));
                
                %acceptance/rejection step
                r=rand(1,1);
                if r<exp(acceptance-logold_IR_dens)*propratio
                    
                    %tracking the acceptances
                    no_of_acceptances(5)=no_of_acceptances(5)+1;
                    
                    %augmented data
                    IR=IR_temp;
                    transtype=transtype_temp;
                    the_hh_with_transition=the_hh_with_transition_temp;
                    
                    %path integrals
                    path_b_infbit=path_b_infbit_temp;
                    path_w_infbit=path_w_infbit_temp;
                    pathrecbit=pathrecbit_temp;
                    timeintervals=timeintervals_temp;
                    
                    %likelihood coefficient terms
                    rmix=rmix_temp;
                    wmix=wmix_temp;
                    bmix=bmix_temp;
                    
                    %states at each transition
                    numinf=numinf_temp;
                    sushh=sushh_temp;
                    infhh=infhh_temp;
                end
            end
        end
    end
    
    %storing samples after burnin
    if jj>burnin
        mcmc_samps(jj-burnin,:)=[bet_new,alph_new,gam_new];
    end
    
    %Prints means of samples, the number of acceptences and the runtime every "update_its" iterations
    if mod(jj,update_its)==0
        jj
        if jj>(burnin+update_its)
            mean(mcmc_samps(jj-burnin-update_its+1:(jj-burnin),:))
            no_of_acceptances
        else
            [bet_new,alph_new,gam_new]
            no_of_acceptances
        end
        toc;
    end
    
    %
end


runtime=toc/3600;

mean(mcmc_samps)
median(mcmc_samps)
cov(mcmc_samps)
end