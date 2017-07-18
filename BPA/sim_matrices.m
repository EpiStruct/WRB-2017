% calcualte the two matrices needed for an efficient stochastic sim of the
% hh model.

function [L,M] = sim_matrices(par)

    N = 3;
   
    beta = par(1);
    gamma = par(2);
    alpha = par(3);
    
    Z = state_map(N);

    % calcualte the l matrix that encodes how each event changes the state of
    % the system.

    % external infection
    l_ext = sparse(2,1,1,10,1);
    
    a2 = (Z(:,1)-Z(:,2));
    a1 = (N-Z(:,1));

    % recovery events
    rows1 = find(a2);
    w = length(rows1);
    l1 = sparse(rows1,1:w,-ones(w,1),10,w);

    states_to = bsxfun(@plus,Z(rows1,:),[0,1]);

    [~,~,rows2] = intersect(states_to,Z,'rows','stable');

    l2 = sparse(rows2,1:w,ones(w,1),10,w);

    % internal infection.
    rows1 = find(a1.*a2>0);
    w = length(rows1);
    l3 = sparse(rows1,1:w,-ones(w,1),10,w);

    states_to = bsxfun(@plus,Z(rows1,:),[1,0]);

    [~,~,rows2] = intersect(states_to,Z,'rows','stable');

    l4 = sparse(rows2,1:w,ones(w,1),10,w);


    l_mat = [l_ext,l3+l4, l1+l2];
    L = [l_mat;sparse([1,1,1,1,-1,-1,-1,-1,-1,-1])];

    % add a one in the H(3,0,0) spot to count total new households.    
    L = L+sparse(1,1,1,11,10);
    
    %/full(l_mat)

    % create the matricies for forming the propensity vector.

    a1a2 = a1.*a2;

    M1 = sparse(2:4,find(a1a2>0),beta*a1a2(a1a2>0)/(N-1),10,11);
    M2 = sparse(1,11,alpha,10,11);
    M3 = sparse(5:10,find(a2>0),gamma*a2(a2>0),10,11);

    M = M1+M2+M3;

end
