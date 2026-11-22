function y = hnf_zero_phase(sos, gain, pad, x)
% Forward-backward filtering in second-order sections. Column boundaries
% are epoch boundaries. Reflect each end once, as in direct-form filtfilt.
% A high-order transfer polynomial can become unstable even for a perfectly
% stable narrowband filter, notably at low EEG frequencies and high Fs.
if size(x,1) <= pad
    error('NeuroFreq:ShortFilterEpoch', ...
        'A %d-sample epoch needs more than %d samples for this filter order.',size(x,1),pad);
end
x = [bsxfun(@minus,2*x(1,:),x(pad+1:-1:2,:)); x; ...
    bsxfun(@minus,2*x(end,:),x(end-1:-1:end-pad,:))];
% Distribute the scalar gain over the cascade so early intermediate sections
% do not amplify a long low-frequency trace before a tiny final gain.
sos(:,1:3) = sos(:,1:3)*abs(gain)^(1/size(sos,1));
sos(1,1:3) = sign(gain)*sos(1,1:3);
for direction = 1:2
    for section = 1:size(sos,1)
        b = sos(section,1:3)/sos(section,4);
        a = sos(section,4:6)/sos(section,4);
        % Constant-input steady-state delay values for a transposed
        % direct-form II biquad, scaled by the first sample of each trace.
        dc = sum(b)/sum(a);
        zi = [dc-b(1); b(3)-a(3)*dc]*x(1,:);
        x = filter(b,a,x,zi,1);
    end
    x = flipud(x);
end
y = x(pad+1:end-pad,:);
end
