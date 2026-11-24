function data = hnf_tf_input(data, Fs)
% A vector represents one sensor; a matrix is sensors x time. Trials always
% occupy dimension 3 and must never share a filter or convolution boundary.
if isa(data,'gpuArray')
    if isempty(data) || issparse(data) || ~any(strcmp(classUnderlying(data), ...
            {'single','double','int8','uint8','int16','uint16','int32','uint32','int64','uint64'}))
        error('NeuroFreq:InvalidTFData','Expected a nonempty dense numeric GPU array.');
    end
else
    validateattributes(data, {'numeric'}, {'nonempty','nonsparse'});
end
validateattributes(Fs, {'numeric'}, {'real','scalar','finite','positive'});
if isvector(data)
    data = reshape(data, 1, []);
end
if ndims(data) > 3 || size(data,2) < 2
    error('NeuroFreq:InvalidTFData', 'Expected sensors x at least two samples x trials.');
end
% FFTs and filters take floating-point inputs; conversion also prevents
% saturation when multiplying integer EEG values by a taper.
if isa(data,'gpuArray')
    floating = any(strcmp(classUnderlying(data),{'single','double'}));
else
    floating = isfloat(data);
end
if ~floating, data = double(data); end
end
