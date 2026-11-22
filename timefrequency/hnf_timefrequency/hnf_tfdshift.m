function out = hnf_tfdshift(in)
% tfdshift -- Shift the spectrum of a TFD by pi radians.
%
%  Usage
%    out = tfdshift(in)
%
%  Inputs
%    in   time-frequency distribution
%
%  Outputs
%    out  shifted time-frequency distribution

% Copyright (C) -- see DiscreteTFDs/Copyright

narginchk(1,1);
out = fftshift(in,1);
end
