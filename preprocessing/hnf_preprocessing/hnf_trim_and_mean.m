function valore=hnf_trim_and_mean(vettore)
% HNF_TRIM_AND_MEAN   Computes average value from vector 'vettore' after removing the top 1% of the values
%
% Usage:
%   >> valore=trim_and_mean(vettore);
%
% Inputs:
%   vettore    - row vector
%
% Outputs:
%   valore     - result
% Count observations independently of vector orientation.
dim = floor(.01*numel(vettore));
tmp=sort(vettore);
valore= mean (tmp(1:(length(vettore)-dim)));


