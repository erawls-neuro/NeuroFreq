function value = hnf_roblerp(left,right,weight)
% HNF_ROBLERP  Linear interpolation of floating-point values, 0 <= weight <= 1.
% Scalar expansion and compatible array sizes follow MATLAB's element-wise
% rules. For finite endpoints each real/imaginary component stays within
% its endpoint range. Exact knots return the supplied endpoint. Intermediate
% values with nonfinite endpoints follow IEEE arithmetic.
if ~isreal(weight) || any(~isfinite(weight(:)) | weight(:)<0 | weight(:)>1)
    error('hnf_roblerp:Weight','Weights must be finite real fractions from zero to one.');
end
if ~isreal(left) || ~isreal(right)
    value = complex(interpolate(real(left),real(right),weight), ...
                    interpolate(imag(left),imag(right),weight));
else
    value = interpolate(left,right,weight);
end
end

function value = interpolate(left,right,weight)
% Start from the nearest endpoint. With opposite signs an endpoint
% difference may overflow even though every convex combination fits.
delta = right-left;
value = left + weight.*delta;
nearRight = weight > .5;
if any(nearRight(:))
    fromRight = right - (1-weight).*delta;
    select = nearRight | false(size(value));
    value(select) = fromRight(select);
end
unsafe = ~isfinite(delta) & isfinite(left) & isfinite(right);
if any(unsafe(:))
    blend = (1-weight).*left + weight.*right;
    select = unsafe | false(size(value));
    value(select) = blend(select);
end
% Returning the operands also avoids 0*Inf at an exact knot.
for side=0:1
    select = (weight==side) | false(size(value));
    if any(select(:))
        if side==0, endpoint=left; else, endpoint=right; end
        endpoint = endpoint + zeros(size(value),'like',value);
        value(select) = endpoint(select);
    end
end
end
