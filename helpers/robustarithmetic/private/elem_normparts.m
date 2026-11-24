function [scale, root, nanSlice, infSlice] = elem_normparts(x,dim)
scale = max(max(abs(real(x)),abs(imag(x))),[],dim);
safe = scale; safe(safe == 0 | ~isfinite(safe)) = 1;
z = x./safe;
root = sqrt(sum(real(z).^2 + imag(z).^2,dim));
nanSlice = any(isnan(x),dim); infSlice = any(isinf(x),dim);
end
