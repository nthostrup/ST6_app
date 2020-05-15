function [res] = pmus2(k,R,C,pmax)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if(k<2.5)
    res = pmax*(1-exp(-(1/(R*C)*k)));
else
    res = pmax*(exp(-(1/(R*C))*(k-2.5)));
end

end

