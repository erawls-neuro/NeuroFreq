function Results = hnf_ept_TFCE_corr_par(data1, data2, chanlocs, varargin)
% Use the identical statistic, spatial adjacency and null model in both modes.
Results = hnf_ept_TFCE_corr(data1, data2, chanlocs, varargin{:}, 'parallel', true);
end
