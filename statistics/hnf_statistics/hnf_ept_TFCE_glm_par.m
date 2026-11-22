function Results = hnf_ept_TFCE_glm_par(data1, data2, chanlocs, varargin)
% Use the identical statistic, spatial adjacency and null model in both modes.
Results = hnf_ept_TFCE_glm(data1, data2, chanlocs, varargin{:}, 'parallel', true);
end
