% Simulates the the full branching process that results from a single hh
% infected at a uniform time in day one.

% H(1) counts number of external infections. 
% H(end) counts the total number of I.

function bh = sim_one_hh(par,REPS)

    % generate stuff needed for the sim.
    [L,M] = sim_matrices(par);
    total_inf = zeros(REPS,1);

    for reps = 1:REPS

        H = zeros(11,1);      % number of hhs in each config.
        H = H + L(:,1);  % inital condition, 1 hh with 1 infected
        t = rand(1);     % start at random time in the first day (0,1]

        while 1==1
            
            a = cumsum(M*H);
            t = t + (1/a(end))*log(1/rand(1));

            if t > 1
                break;
            end

            index = find(a>=rand(1)*a(end),1,'first');
            H = H + L(:,index);
            a = cumsum(M*H);
            t = t + (1/a(end))*log(1/rand(1));

        end

        total_inf(reps) = H(1);
    end

    bh = hist(total_inf,0:50)/REPS;
    
end