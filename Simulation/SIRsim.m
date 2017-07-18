%Simulates SIR stochastic household epidemics

% household size
k=3;
% number of households
m=50000;

%initilising a while loop to ensure the epidemic does not die out
counter=1;
infected=0;
while infected==0 
    
    % initial condition
    delWmat=zeros(hthingo+100,200);
    delWmat(1,2)=1;
    Initialstate=zeros(m,2); %S I
    Initialstate(2:m,1)=k;
    Initialstate(1,1)=k-1;
    Initialstate(1,2)=1;
    
    % pre-specifying vector sizes for efficiency
    BIGNUM=10000;
    time=[rand(1,1),zeros(1,BIGNUM-1)];%uniform on first day
    infected=1;
    totalhh=[1,zeros(1,BIGNUM-1)];
    beeninfected=[1,zeros(1,m-1)];
    delHt=[0,1,zeros(1,198)];
    ordering=zeros(1,hthingo+100);
    ordering(1)=1;
    
    counter=1;
    counter2=1;
    Currstate=Initialstate;
    r=zeros(m,2);
    
    while infected>0
        
        %rates of each possible transition
        tiff=(sum(Currstate(:,2))-Currstate(:,2))/(k*(m-1));
        r(:,1)=bet*(Currstate(:,2).*Currstate(:,1))./(k-1)+alph*(tiff.*Currstate(:,1)); %infection event rate in each household
        r(:,2)=gam*(Currstate(:,2)); %recovery rate in each household
        
        %Exponentially distributed random time until the next transition
        leave=exprnd(1/sum(sum(r))); 
        
        %Generating the transition type and household
        whichtransition=r(:)./sum(sum(r));
        cumulative=cumsum(whichtransition);
        random=rand(1,1);
        j=sum((cumulative-random)<0)+1; %which transition specified
        l=0:1;
        b=sum((j-l*m)>0); %the transition type
        a=j-m*(b-1);   %the household
        
        %Updating the state
        Currstate(a,b)=Currstate(a,b)-1;
        if b==1
            Currstate(a,b+1)=Currstate(a,b+1)+1;
        end
        counter=counter+1;
        
        %Keeping track of the exact times of events and the corresponding
        %number of infectives and cumulative number of infection events
        time(counter)=time(counter-1)+leave;
        infected=sum(Currstate(:,2));
        
        %Saves the time at which the specified household level of infection is
        %reached and breaks out of the while loop after the day has finished
        %being simulated
        if totalhh(counter-1)>=hthingo
            timetostop=ceil(time(totalhh==hthingo));
            T_horizon=timetostop(1);
            if time(counter)>=T_horizon
                break
            end
        end
        
        %Saves FF100 type data 
        if b==1
            
            %If a new household is infected
            if beeninfected(a)==0
                %keep track of the number of infected households
                counter2=counter2+1;
                
                %a vector that keeps track of the order in which households
                %became infected
                ordering(counter2)=a;
                
                %matrix of the number of infections in each household on
                %each day
                delWmat(counter2,ceil(time(counter))+1)=1;
                
                %vector of the total number of newly infected households on
                %each day
                delHt(ceil(time(counter))+1)=delHt(ceil(time(counter))+1)+1;
                
                %keeping track of the number of households at each
                %transition time
                totalhh(counter)=totalhh(counter-1)+1;
                
                %noting that household "a" has now been infected
                beeninfected(a)=1;
            else
                %if the household has already been infected save the
                %relevant infections in the relevant section of the data
                %matrix
                delWmat(ordering==a,ceil(time(counter))+1)=delWmat(ordering==a,ceil(time(counter))+1)+1;
                totalhh(counter)=totalhh(counter-1);
            end
        else
            totalhh(counter)=totalhh(counter-1);
        end
    end
end

%Keeping the relevant bits of data matricies
Wmat=cumsum(delWmat(1:counter2,1:(T_horizon+1)),2);
delHt=delHt(1:(T_horizon+1));
Ht=zeros(1,T_horizon+1);
for t=1:(T_horizon+1)
   Ht(t)=sum(delHt(1:t)); 
end

Wmat=Wmat(:,2:end);
Ht=Ht(2:end);
