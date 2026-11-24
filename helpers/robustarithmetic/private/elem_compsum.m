function total = elem_compsum(x,dim,mask)
% Internal: scaled Neumaier accumulation over selected finite slices.
if dim > ndims(x), total = x; return; end
order = [dim,1:dim-1,dim+1:ndims(x)];
v = permute(x,order); v = reshape(v,size(x,dim),[]); v = v(:,mask(:));
total = elem_realcompsum(real(v));
if ~isreal(x), total = complex(total,elem_realcompsum(imag(v))); end
end

function y = elem_realcompsum(v)
s = max(abs(v),[],1); s(s == 0) = 1;
% Scaling can discard small but representable observations after two large
% contributions cancel. Keep all significand bits for those rare columns.
magnitudes = abs(v); magnitudes(magnitudes==0) = Inf;
exact = min(magnitudes,[],1)./s < realmin(class(v))/eps(class(v));
if any(exact)
    y = zeros(1,size(v,2),'like',v);
    for col=find(exact), y(col) = elem_exactsum(v(:,col)); end
    if any(~exact), y(~exact) = elem_realcompsum(v(:,~exact)); end
    return
end
if isa(v,'double') && numel(v)>=4096 && ...
        exist('hnf_mex_available','file')==2 && hnf_mex_available('hnf_compsum_mex')
    try
        y = hnf_compsum_mex(v,s);
        return
    catch err
        if ~hnf_mex_load_error(err), rethrow(err); end
    end
end
v = v./s;
t = zeros(1,size(v,2),'like',v); c = t;
for k = 1:size(v,1)
    a = v(k,:); u = t + a;
    choose = abs(t) >= abs(a);
    delta = (a-u)+t; delta(choose) = (t(choose)-u(choose))+a(choose);
    c = c + delta; t = u;
end
y = (t+c).*s;
end
