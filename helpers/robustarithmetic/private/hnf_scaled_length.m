function [y,info] = hnf_scaled_length(x,dim,nanflag,rmsMode)
% Octave/MATLAB NaN conventions for ordinary Euclidean norm and RMS.
% The general policy forms of hnf_robnorm and hnf_robrms remain separate.
elem_float(x,'hnf_scaled_length');
dim=hnf_aux_dim(x,dim);
if ~any(strcmp(nanflag,{'omitnan','includenan'}))
    error('hnf:Options','Use includenan or omitnan.');
end
omit=strcmp(nanflag,'omitnan');
if nargout==1 && isreal(x) && numel(x)>=65536 && ...
        exist('hnf_mex_available','file')==2 && hnf_mex_available('hnf_scaled_length_mex')
    try
        y=hnf_scaled_length_mex(x,double(dim),omit,rmsMode);
        return
    catch problem
        if ~hnf_mex_load_error(problem), rethrow(problem); end
    end
end
nanMask=isnan(x);
infSlice=any(isinf(x) & ~nanMask,dim);
nanSlice=any(nanMask,dim) & ~omit;
n=size(x,dim);
if omit, count=n-sum(nanMask,dim); else, count=n; end
x(~isfinite(x))=0;
if n==0
    y=sum(real(x),dim);
    if rmsMode, y(:)=NaN; end
else
    scale=max(abs(real(x)),[],dim);
    if ~isreal(x), scale=max(scale,max(abs(imag(x)),[],dim)); end
    scale(scale==0)=1;
    z=x./scale;
    squares=real(z).^2;
    if ~isreal(z), squares=squares+imag(z).^2; end
    q=sum(squares,dim,'double');
    if rmsMode, q=q./max(count,1); end
    y=cast(double(scale).*sqrt(q),'like',real(x));
    if rmsMode && omit, y(count==0)=NaN; end
end
y(infSlice)=Inf; y(nanSlice)=NaN;
if nargout>1
    info=struct('dimension',dim,'count',count,'nonfiniteOutput',nnz(~isfinite(y)), ...
        'nanPolicy',nanflag);
end
end
