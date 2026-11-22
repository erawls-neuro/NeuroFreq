function data = hnf_tf_input(data, Fs)
% A vector represents one sensor; a matrix is sensors x time. Trials always
% occupy dimension 3 and must never share a filter or convolution boundary.
validateattributes(data, {'numeric'}, {'nonempty','nonsparse'});
validateattributes(Fs, {'numeric'}, {'real','scalar','finite','positive'});
if isvector(data)
    data = reshape(data, 1, []);
end
if ndims(data) > 3 || size(data,2) < 2
    error('NeuroFreq:InvalidTFData', 'Expected sensors x at least two samples x trials.');
end
% FFTs and filters take floating-point inputs; conversion also prevents
% saturation when multiplying integer EEG values by a taper.
if ~isfloat(data), data = double(data); end
end
