function [a,b] = elem_pair(a,b,id)
elem_float(a,id); elem_float(b,id);
if strcmp(class(a),class(b)), return; end
if isscalar(a) && ~isscalar(b)
    before = a; a = cast(a,'like',b);
    elem_castcheck(before,a,id);
elseif isscalar(b) && ~isscalar(a)
    before = b; b = cast(b,'like',a);
    elem_castcheck(before,b,id);
else
    error([id ':Type'],'Use matching classes, or one scalar of the other floating-point class.');
end
end

function elem_castcheck(before,after,id)
if (isfinite(before) && ~isfinite(after)) || ...
        (real(before) ~= 0 && real(after) == 0) || ...
        (imag(before) ~= 0 && imag(after) == 0)
    error([id ':Type'],'Mixed-class scalar conversion would overflow or underflow; cast explicitly.');
end
end
