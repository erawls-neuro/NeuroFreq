function guess = hnf_specparam_slope_guess(freqs,spectrum)
% Endpoint slope in log10-Hz. Take a log difference for very wide grids
% and log1p for adjacent frequencies whose separately rounded logs coincide.
lo = double(freqs(1));
hi = double(freqs(end));
if hi-lo <= lo
    decades = log1p((hi-lo)/lo)/log(10);
else
    decades = log10(hi)-log10(lo);
end
first = double(spectrum(1)); last = double(spectrum(end));
change = first-last;
if isinf(change) && isfinite(first) && isfinite(last)
    guess = (first/2-last/2)/(decades/2);
else
    guess = change/decades;
end
end
