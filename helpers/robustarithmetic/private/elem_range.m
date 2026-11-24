function y = elem_range(y, finiteSource, nonzeroSource, o, id)
% Internal: policies affect finite-source overflow, never a missing value.
if ~strcmp(o.overflow,'allow')
    over = finiteSource & isinf(y);
    if any(over(:))
        if strcmp(o.overflow,'error')
            error([id ':Overflow'],'The result exceeds the output floating-point range.');
        end
        r = real(y); i = imag(y); limit = realmax(class(y));
        hit = over & isinf(r); r(hit) = sign(r(hit)).*limit;
        if isreal(y)
            y = r;
        else
            hit = over & isinf(i); i(hit) = sign(i(hit)).*limit;
            y = complex(r,i);
        end
    end
end
if strcmp(o.underflow,'allow'), return; end
under = finiteSource & nonzeroSource & (y == 0);
if any(under(:))
    if strcmp(o.underflow,'error')
        error([id ':Underflow'],'A nonzero result rounded to zero.');
    end
    % The caller supplies the mathematically meaningful sign/phase when
    % needed. This floor is only used directly for nonnegative results.
    y(under) = realmin(class(y)).*eps(class(y));
end
end
