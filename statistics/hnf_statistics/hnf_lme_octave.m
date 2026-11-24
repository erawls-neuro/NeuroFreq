function [beta,tstat] = hnf_lme_octave(Y,X,subject,slopes,method,diagonal,backend)
% Fit the numeric Gaussian mixed model used by hnf_ept_TFCE_lme.
% Y is observations x locations. Return non-intercept fixed effects in
% design-column order. A diagonal random covariance requires a separate
% independent variance component for each column, including the intercept.
fixed = [ones(size(X,1),1), X];
if rank(fixed) < size(fixed,2)
    error('NeuroFreq:LMEDesign','Fixed-effect design with intercept is rank deficient.');
end
if ~any(strcmpi(method,{'ML','REML'}))
    error('NeuroFreq:LMEFitMethod','Choose ML or REML.');
end
if strcmp(backend,'R lme4 lmer')
    [beta,tstat] = local_fit_lme4(Y,X,subject,slopes,method,diagonal);
    return
end
Z = [fixed(:,1),X(:,slopes)];
group = subject(:);
if diagonal && size(Z,2) > 1
    Z = num2cell(Z,1);
    group = repmat({group},size(Z));
end
beta = zeros(size(X,2),size(Y,2));
tstat = beta;
for k = 1:size(Y,2)
    model = fitlmematrix(fixed,Y(:,k),Z,group,'FitMethod',method);
    coefficients = model.Coefficients;
    beta(:,k) = coefficients.Estimate(2:end);
    tstat(:,k) = coefficients.tStat(2:end);
end
end

function [beta,tstat] = local_fit_lme4(Y,X,subject,slopes,method,diagonal)
% One R process per design (all voxels), with full-precision binary I/O.
% Fixed command text and environment variables keep paths out of shell syntax.
folder = tempname;
[ok,msg] = mkdir(folder);
if ~ok, error('NeuroFreq:LMEIO','%s',msg); end
cleanFolder = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
fid = fopen(fullfile(folder,'data.bin'),'wb','ieee-le');
if fid < 0, error('NeuroFreq:LMEIO','Cannot write the temporary model input.'); end
closeInput = onCleanup(@() fclose(fid));
header = [size(Y,1), size(Y,2), size(X,2), numel(slopes), strcmpi(method,'REML'), diagonal];
fwrite(fid,header,'int32');
fwrite(fid,slopes,'int32');
fwrite(fid,subject,'double');
fwrite(fid,X,'double');
fwrite(fid,Y,'double');
clear closeInput

oldFolder = getenv('NEUROFREQ_LME_DIR');
oldScript = getenv('NEUROFREQ_LME_SCRIPT');
envCleanup = onCleanup(@() local_restore_env(oldFolder,oldScript)); %#ok<NASGU>
setenv('NEUROFREQ_LME_DIR',folder);
setenv('NEUROFREQ_LME_SCRIPT',fullfile(fileparts(mfilename('fullpath')),'hnf_lme_lme4.R'));
[status,output] = system('Rscript --vanilla -e "source(Sys.getenv(''NEUROFREQ_LME_SCRIPT''))" 2>&1');
if status ~= 0
    error('NeuroFreq:LMEBackend',['Octave LME requires statistics fitlmematrix or Rscript with lme4. ' ...
        'R backend failed:\n%s'],output);
end
if ~isempty(strtrim(output)), warning('NeuroFreq:LMEBackendWarning','%s',strtrim(output)); end
fid = fopen(fullfile(folder,'fit.bin'),'rb','ieee-le');
if fid < 0, error('NeuroFreq:LMEIO','R returned without a coefficient file.'); end
closeOutput = onCleanup(@() fclose(fid)); %#ok<NASGU>
[both,n] = fread(fid,[size(X,2),2*size(Y,2)],'double');
if n ~= 2*size(X,2)*size(Y,2)
    error('NeuroFreq:LMEIO','R returned an incomplete coefficient matrix.');
end
beta = both(:,1:size(Y,2));
tstat = both(:,size(Y,2)+1:end);
end

function local_restore_env(folder,script)
setenv('NEUROFREQ_LME_DIR',folder);
setenv('NEUROFREQ_LME_SCRIPT',script);
end
