function valore=hnf_trim_and_max(vettore)
% HNF_TRIM_AND_MAX   Computes max value from vector 'vettore' after removing the top 1% of the values
%
% Usage:
%   >> valore=trim_and_max(vettore);
%
% Inputs:
%   vettore    - row vector
%
% Outputs:
%   valore     - result
% Count observations independently of vector orientation.
dim = floor(.01*numel(vettore));
tmp=sort(vettore);
valore= tmp(length(vettore)-dim);


