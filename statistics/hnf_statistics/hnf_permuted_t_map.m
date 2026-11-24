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
    signed = hnf_stats_scale(data) .* draw;
    origin = signed(1,:,:,:);
    centered = signed-origin;
    t = hnf_robmean(signed,1,'includenan') ./ (std(centered,0,1) / sqrt(size(data,1)));
else
    data = hnf_stats_scale(data);
    a = data(draw(1:nA),:,:,:);
    b = data(draw(nA+1:end),:,:,:);
    origA = a(1,:,:,:); origB = b(1,:,:,:);
    [a,scaleA] = hnf_stats_scale(a-origA);
    [b,scaleB] = hnf_stats_scale(b-origB);
    diffMean = (origA-origB) + (mean(a,1).*scaleA - mean(b,1).*scaleB);
    seA = (std(a,0,1)/sqrt(nA)) .* scaleA;
    seB = (std(b,0,1)/sqrt(size(b,1))) .* scaleB;
    t = diffMean ./ hypot(seA,seB);
end
end
