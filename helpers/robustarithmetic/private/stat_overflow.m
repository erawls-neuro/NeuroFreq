function [y, overflow] = stat_overflow(y, policy)
overflow = isinf(real(y)) | isinf(imag(y));
if any(overflow(:))
    switch policy
        case 'error'
            error('hnf:statistics:Overflow','The requested statistic is not representable in the input class.');
        case 'saturate'
            if isreal(y)
                y(isinf(y)) = sign(y(isinf(y)))*realmax(class(y));
            else
                re = real(y); im = imag(y);
                re(isinf(re)) = sign(re(isinf(re)))*realmax(class(y));
                im(isinf(im)) = sign(im(isinf(im)))*realmax(class(y));
                y = complex(re,im);
            end
    end
end
end
