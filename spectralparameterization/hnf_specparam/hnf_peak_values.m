function [values,basis,u] = hnf_peak_values(freqs,centers,heights,widths,shape)
% Evaluate rows of Gaussian or Lorentzian peaks, including their tails.
distance = bsxfun(@minus,freqs,centers);
u = bsxfun(@rdivide,distance,widths);
if strcmp(shape,'cauchy')
    basis = 1./(1+u.^2);
    values = bsxfun(@times,heights,basis);
    % Squaring u can overflow while height/u² remains representable.
    tail = abs(u) > sqrt(realmax)/2;
    if any(tail(:))
        [rows,~] = find(tail);
        d = distance(tail);
        w = widths(rows);
        h = heights(rows);
        inv = w(:)./d(:);
        basis(tail) = inv.^2./(1+inv.^2);
        values(tail) = ((h(:).*inv).*inv)./(1+inv.^2);
    end
else
    e = -0.5.*u.*u;
    basis = exp(e);
    values = bsxfun(@times,heights,basis);
    % A very large height can offset underflow of the unweighted tail.
    tail = bsxfun(@and,e < log(realmin),abs(heights)>1);
    if any(tail(:))
        [rows,~] = find(tail);
        h = heights(rows);
        tailLog = e(tail);
        values(tail) = sign(h(:)).*exp(log(abs(h(:)))+tailLog(:));
    end
end
end
