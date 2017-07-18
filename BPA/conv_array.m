% A function for calculating the convulution matrix
% Inputs:
% day_dist- the pmf of the size of a household branching process over a day 
% number- the number of columns in the convolution matrix
% cutoff- the number of rows in the convolution matrix
%
% Output:
% M- the convolution matrix
function M = conv_array(day_dist,number,cutoff)

    M = zeros(cutoff,number);
    
    % for zero infections
    M(1,1) = 1; 
    
    ci = day_dist';
    M(1:length(ci),2) = ci;

    for ii=3:number
        ci = conv(ci,day_dist');

        if length(ci) >= cutoff
            M(:,ii) = ci(1:cutoff);
        else
            M(1:length(ci),ii) = ci;
        end

    end

end