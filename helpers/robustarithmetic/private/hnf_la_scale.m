function s = hnf_la_scale(A)
% Internal: finite component scale (complex magnitude itself may overflow).
if isempty(A), s = cast(1,'like',real(A)); return; end
s = max(abs(real(A(:))));
if ~isreal(A), s = max(s,max(abs(imag(A(:))))); end
if s == 0, s = cast(1,'like',real(A)); end
end
