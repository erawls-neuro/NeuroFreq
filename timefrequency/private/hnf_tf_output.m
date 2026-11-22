function output = hnf_tf_output(values)
% Keep frequency, time, and trial axes intact, even when only one frequency
% or one STFT frame is requested. NF omits only a singleton sensor axis.
if size(values,1) == 1
    output = reshape(values, size(values,2), size(values,3), size(values,4));
else
    output = values;
end
end
