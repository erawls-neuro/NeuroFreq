function x = elem_input(x, o, id)
if strcmp(o.inputPolicy,'propagate'), return; end
bad = ~isfinite(x);
if ~any(bad(:)), return; end
if strcmp(o.inputPolicy,'error')
    error([id ':NonfiniteInput'],'Input contains NaN or Inf.');
end
fill = cast(o.fillValue,'like',x);
if ~isfinite(fill)
    error([id ':Options'],'fillValue is not representable in the input class.');
end
x(bad) = fill;
end
