function ys = hnf_gaussian(freqs, mu, hgt, sigma)
%       Gaussian function to use for fitting.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency vector to create gaussian fit for.
%       mu, hgt, sigma : doubles
%           Parameters that define gaussian function (centre frequency,
%           height, and standard deviation).
%
%       Returns
%       -------
%       ys :    1xn array
%       Output values for gaussian function.

    if isscalar(mu) && isscalar(hgt) && isscalar(sigma)
        ys = reshape(hnf_peak_values(freqs(:).',mu,hgt,sigma,'gaussian'),size(freqs));
    else
        u = (freqs-mu)./sigma;
        ys = hgt.*exp(-.5*u.^2);
    end

end
