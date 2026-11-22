function t = hnf_permuted_t_map(data, draw, useNative, nA)
% The outer permutation loop draws the random assignment once. The MEX
% and MATLAB backends receive the same assignment, including on fallback.
% Omitting nA selects a sign flip. Supplying it selects the two-group test.
% Serial loops pass the local eligibility decision. Parallel loops pass
% the exact MEX file selected on the client (empty when native is off).
if ischar(useNative)
    useNative = ~isempty(useNative) && exist('hnf_permuted_t_mex','file') == 3 && ...
        strcmp(which('hnf_permuted_t_mex'),useNative);
end
if useNative
    try
        if nargin == 3
            t = hnf_permuted_t_mex(data,draw);
        else
            t = hnf_permuted_t_mex(data,draw,nA);
        end
        return
    catch failure
        % The client may see the binary while a parallel worker does not.
        % A missing worker binary can use the already-drawn assignment.
        unavailable = exist('hnf_permuted_t_mex','file') ~= 3;
        if ~unavailable && ~hnf_mex_load_error(failure), rethrow(failure); end
        persistent warned
        if isempty(warned)
            warning('NeuroFreq:NativeTLoad','%s Using MATLAB for permutation t statistics.',failure.message);
            warned = true;
        end
    end
end
if nargin == 3
    signed = draw .* data;
    t = mean(signed,1) ./ (std(signed,0,1) ./ sqrt(size(data,1)));
else
    a = data(draw(1:nA),:,:,:);
    b = data(draw(nA+1:end),:,:,:);
    t = (mean(a)-mean(b)) ./ sqrt(var(a)/nA + var(b)/size(b,1));
end
end
