function [peaks, locations, widths, prominence] = hnf_mara_peaks(y,x)
% MARA uses signed, finite spectra, sorted coordinates and half-prominence widths.
if ~exist('OCTAVE_VERSION','builtin')
    [peaks,locations,widths,prominence] = findpeaks(y,x);
    return
end
% Use the first sample of a flat maximum. Endpoints cannot be peaks.
starts = [1 find(diff(y)~=0)+1];
levels = y(starts);
at = find(levels(2:end-1)>levels(1:end-2) & levels(2:end-1)>levels(3:end))+1;
idx = starts(at);
peaks = y(idx);
locations = x(idx);
widths = zeros(size(peaks));
prominence = zeros(size(peaks));
for k = 1:numel(idx)
    i = idx(k);
    % Extend horizontally until a higher sample or the edge on each side.
    left = find(y(1:i)>y(i),1,'last');
    right = find(y(i:end)>y(i),1)+i-1;
    if isempty(left), left = 1; end
    if isempty(right), right = numel(y); end
    floorLevel = max(min(y(left:i)),min(y(i:right)));
    prominence(k) = y(i)-floorLevel;
    half = (y(i)+floorLevel)/2;
    a = find(y(1:i)<=half,1,'last');
    b = find(y(i:end)<=half,1)+i-1;
    xLeft = x(a)+(half-y(a))/(y(a+1)-y(a))*(x(a+1)-x(a));
    xRight = x(b)+(half-y(b))/(y(b-1)-y(b))*(x(b-1)-x(b));
    widths(k) = xRight-xLeft;
end
end
