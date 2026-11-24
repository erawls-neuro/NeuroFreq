function [s, mu] = stat_rawmad(raw, valid)
% Order raw data before arithmetic scaling: remote outliers must not erase MAD.
mu = stat_quantile(raw,valid,.5);
residual = abs(bsxfun(@minus,raw,mu));
s = stat_quantile(residual,valid,.5);
end
