function z = elem_cdivfinite(a,b)
% Finite, nonzero-denominator complex division by component exponents.
% Cross products are added at a common exponent before final scaling;
% tiny numerator components are not discarded by whole-number scaling.
[ar,ear] = log2(abs(real(a))); ar = ar.*sign(real(a));
[ai,eai] = log2(abs(imag(a))); ai = ai.*sign(imag(a));
[br,ebr] = log2(abs(real(b))); br = br.*sign(real(b));
[bi,ebi] = log2(abs(imag(b))); bi = bi.*sign(imag(b));
er = ebr; ei = ebi; er(br == 0) = -Inf; ei(bi == 0) = -Inf;
eb = max(er,ei);
ebr(br == 0) = eb(br == 0); ebi(bi == 0) = eb(bi == 0);
dr = pow2(br,ebr-eb); di = pow2(bi,ebi-eb);
d = dr.^2 + di.^2;
[nr,enr] = elem_addproducts(ar.*br,ear+ebr,ai.*bi,eai+ebi);
[ni,eni] = elem_addproducts(ai.*br,eai+ebr,-ar.*bi,ear+ebi);
enr = enr-2.*eb; eni = eni-2.*eb;
enr(nr == 0) = 0; eni(ni == 0) = 0;
z = complex(pow2(nr./d,enr),pow2(ni./d,eni));
end

function [n,e] = elem_addproducts(p,ep,q,eq)
zp = p == 0; zq = q == 0;
ep(zp) = -Inf; eq(zq) = -Inf;
e = max(ep,eq); e(zp & zq) = 0;
% Give zero terms a benign exponent instead of passing -Inf to pow2.
ep(zp) = e(zp); eq(zq) = e(zq);
n = pow2(p,ep-e) + pow2(q,eq-e);
end
