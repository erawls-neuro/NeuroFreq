function elem_float(x, id)
% Internal: explicitly constrain the supported storage and arithmetic types.
if ~(strcmp(class(x),'double') || strcmp(class(x),'single')) || issparse(x)
    error([id ':Type'],'Inputs must be full CPU single or double arrays.');
end
end
