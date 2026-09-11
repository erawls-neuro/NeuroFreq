function TFCE_Obs = hnf_tfce_transformation(T_Obs,ChN,E_H,flag_ft)
%
% Transforms a 1/2/3D tensor of rscores to TFCE via enhancing clusters of
% tests. Should only be called internally from the ept_TFCE functions.
%

% check for non-zero T-values (will crash TFCE)
if max(abs(T_Obs(:))) < 0.00001
    error('T-values were all 0')
end
% TFCE transformation...
if ismatrix(T_Obs)
    if ~flag_ft
        TFCE_Obs = ept_mex_TFCE2D(T_Obs, ChN, E_H);
    else
        % if data is not in channel neighbourhood
        % artificially create 3rd dimension
        T_Obs = repmat(T_Obs, [1, 1, 2]);
        deltaT = max(abs(T_Obs(:)))/50;
        TFCE_Obs = ept_mex_TFCE(T_Obs, deltaT);
        % remove extra dimension
        T_Obs = T_Obs(:, :, 1);
        TFCE_Obs = TFCE_Obs(:, :, 1);
    end
end
if ndims(T_Obs) == 3
    TFCE_Obs = ept_mex_TFCE3D(T_Obs, ChN, E_H);
end
end