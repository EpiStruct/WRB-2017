% Calcualte the map between states of the system and event counts.
% SIR model, N=size of system.

function X = state_map(N)
    
    K = (N+1)*(N+2)/2;
    X = zeros(K,2);
    
    count = 1;
    for ii = 0:N
        for jj = 0:ii
            X(count,:) = [ii;jj];
            count = count +1;
        end
    end
    
end