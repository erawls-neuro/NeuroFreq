function [eh, count, grid, neighbours, parallel, corrMode] = hnf_tfce_options(data, chanlocs, varargin)
% Shared option contract for correlation and GLM maps.
p = inputParser;
addParameter(p, 'e_h', [0.66 2]);
addParameter(p, 'nperm', 1000, @(x) isscalar(x) && isfinite(x) && x >= 1 && x == fix(x));
addParameter(p, 'flag_ft', isempty(chanlocs));
addParameter(p, 'chn', []);
addParameter(p, 'parallel', false);
addParameter(p,'corr_mode','auto',@(x) any(strcmpi(x,{'auto','pointwise','predictors'})));
parse(p, varargin{:});
eh = p.Results.e_h;
count = p.Results.nperm;
grid = p.Results.flag_ft;
neighbours = p.Results.chn;
parallel = p.Results.parallel;
corrMode = p.Results.corr_mode;
if ~grid
    if numel(chanlocs) ~= size(data,2)
        error('NeuroFreq:ChanlocMismatch', 'Channel locations must match the first map dimension.');
    end
    if isempty(neighbours)
        if numel(chanlocs) == 1, neighbours = 1; else, neighbours = ept_ChN2(chanlocs); end
    end
end
end
