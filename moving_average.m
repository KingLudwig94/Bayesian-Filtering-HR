function y = moving_average(x,k)
%MOVING_AVERAGE Summary of this function goes here
%   Detailed explanation goes here
    L = length(x);
    y = x;
    for i = ceil(k/2):L-ceil(k/2)
        y(i) = mean(x(i-floor(k/2):i+floor(k/2)));
    end
    
end

