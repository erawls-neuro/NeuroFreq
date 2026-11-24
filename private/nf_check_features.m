function rows = nf_check_features(selectedPaths, root, eeglabRoot)
% Check optional capabilities by feature. MATLAB toolboxes follow the user's
% MATLAB license. They are never downloaded into NF or checked out here.
% Temporarily expose selected providers so the MEX and classifier checks see
% the same files that nf_run will use. No scientific computation is run.
originalPath = path;
restore = onCleanup(@() path(originalPath)); %#ok<NASGU>
if nargin < 2, root = fileparts(fileparts(mfilename('fullpath'))); end
if nargin < 3, eeglabRoot = ''; end
octave = exist('OCTAVE_VERSION','builtin') ~= 0;
if octave, packages = nf_octave_packages(); end
core = nf_eeglab_paths(eeglabRoot);
providers = [nf_source_dirs(root), selectedPaths(:)', core];
addpath(providers{:}, '-begin');

% Each row gives a product, its MATLAB license feature, representative
% entry points, and the NF methods that use it. A license test alone cannot
% detect a licensed product whose files are missing from this installation.
products = {
    'Signal Processing Toolbox', 'Signal_Toolbox', {'butter','filtfilt','hilbert','mscohere','pwelch','findpeaks'}, ...
        'Filtering, analytic signals, and preprocessing quality metrics'
    'Wavelet Toolbox', 'Wavelet_Toolbox', {'cwt','wdenoise'}, ...
        'nf_cwt and HAPPE+ER wavelet precleaning'
    'Statistics and Machine Learning Toolbox', 'Statistics_Toolbox', {'fitlme','robustfit','pdist','classify'}, ...
        'Mixed-effects models, robust regression, and ICA classifiers'
    'Optimization Toolbox', 'Optimization_Toolbox', {'fmincon','lsqcurvefit'}, ...
        'Constrained peak fitting and one of the MARA fit options'
    'Curve Fitting Toolbox', 'Curve_Fitting_Toolbox', {'fit','fitoptions','fittype'}, ...
        'Alternative MARA fitter when lsqcurvefit is unavailable'
    'Parallel Computing Toolbox', 'Distrib_Computing_Toolbox', {'gpuDeviceCount','parpool'}, ...
        'Parallel estimation and GEDAI v1.7 (including its CPU mode)'
    };
if octave
    % Probe the entry points that Octave implementations actually use.
    % An installed package alone does not establish support for a call signature.
    products = {
        'Signal processing', '', {'butter','filtfilt','hilbert','mscohere','pwelch','findpeaks'}, ...
            'Signal package routines used by preprocessing and analysis'
        'Optimization', '', {'fmincon','lsqcurvefit'}, 'Constrained optimization'
        'Table data', '', {'table','categorical','struct2table'}, 'Table operations from datatypes or Tablicious'
        'Octave viewer', '', {'nf_viewerapp'}, 'Figure-based viewer for nf_run(''nf_viewerapp'',...)'
        };
end
rows = cell(size(products,1),3);
for i = 1:size(products,1)
    probes = products{i,3};
    absent = probes(cellfun(@(name) isempty(which(name)), probes));
    if ~isempty(absent)
        status = ['Missing functions: ' strjoin(absent, ', ')];
    elseif octave
        status = 'Entry points found.';
    elseif ~license('test',products{i,2})
        status = 'Functions found; the MATLAB license test did not pass.';
    else
        status = 'Functions found and license test passed.';
    end
    rows(i,:) = {products{i,1},status,products{i,4}};
end

if octave, rows = [packages; rows]; end

% A shared signed integrator serves both runtimes. EPT still supplies the
% optional automatic neighbour graph for electrode montages.
probes = {'hnf_tfce_integrate'};
rows(end+1,:) = {'TFCE enhancement',nf_probe_result(probes(cellfun(@(n) isempty(which(n)),probes))), ...
    'Native or interpreted integration; automatic multi-channel neighbours require ept_ChN2'};
if ~isempty(which('pop_iclabel'))
    icRoot = fileparts(which('pop_iclabel'));
    binary = fullfile(icRoot,'matconvnet','matlab','mex',['vl_nnconv.' mexext]);
    if isfile(binary) || exist('vl_nnconv','file') == 3
        status = ['Found vl_nnconv.' mexext '.'];
    else
        status = sprintf('No vl_nnconv.%s found. A MatConvNet build for this MATLAB architecture may be required.', mexext);
    end
    rows(end+1,:) = {'ICLabel native convolution',status,'ICLabel component classification'};
end
if ~isempty(which('runamica15'))
    if ispc
        binary = 'amica15mkl.exe';
    elseif ismac
        binary = 'amica15mac';
    else
        binary = 'amica15ex';
    end
    if isfile(fullfile(fileparts(which('runamica15')),binary))
        status = ['Found ' binary ' (execution not tested).'];
    else
        status = ['Missing ' binary ' beside runamica15.m.'];
    end
    rows(end+1,:) = {'AMICA executable',status,'AMICA decomposition'};
end

% A distribution can match GitHub perfectly while a particular method still
% has unresolved calls. Report those separately from missing vendor packages.
if ~isempty(eeglabRoot)
    needed = {'eeg_checkset','pop_select','pop_resample','pop_saveset', ...
        'convertlocs','eeg_interp','pop_reref','pop_runica','pop_subcomp', ...
        'pop_eegthresh','pop_rejspec','pop_rejepoch','pop_epoch','eeg_regepochs','pop_rmbase'};
    missing = needed(~cellfun(@(n) nf_under(which(n),eeglabRoot), needed));
    rows(end+1,:) = {'Selected EEGLAB core', nf_probe_result(missing), ...
        'Core tools called by nf_preprocess; individual presets use a subset'};
end

source = nf_source_dirs(root);
rows = [rows; nf_helper_checks(source(2:end))];
classifier = fullfile(root,'preprocessing','nf_cleanic.m');
if isfile(classifier) && ~isempty(regexp(fileread(classifier),'\<nfADJUST\s*\(', 'once'))
    needed = {'nfADJUST','nf_adjust_core','nf_adjust_role_atlas'};
    missing = needed(~cellfun(@(n) nf_under(which(n),root), needed));
    rows(end+1,:) = {'NeuroFreq native ADJUST', nf_probe_result(missing), ...
        'Native nfadjust classification in nf_cleanic'};
end

for i = 1:size(rows,1)
    fprintf('  %s: %s\n    %s\n',rows{i,1},rows{i,2},rows{i,3});
end
end

function status = nf_probe_result(missing)
if isempty(missing)
    status = 'Required functions found.';
else
    status = ['Missing functions: ' strjoin(missing, ', ')];
end
end

function rows = nf_helper_checks(folders)
% Older callers sometimes use an unprefixed helper name, e.g. binomial2
% when this tree ships hnf_binomial2.m. Build the candidates from disk, so
% correcting a call or supplying a compatible function clears the finding.
rows = cell(0,3);
aliases = {};
for k = 1:numel(folders)
    helpers = dir(fullfile(folders{k},'hnf_*.m'));
    aliases = [aliases regexprep({helpers.name},'^hnf_|\.m$','')]; %#ok<AGROW>
end
aliases = aliases(cellfun(@(name) isempty(which(name)), aliases));
if isempty(aliases), return; end
pattern = ['(?<![.\w])(?:' strjoin(aliases,'|') ')(?=\s*\()'];
for k = 1:numel(folders)
    % Keep an entry per directory so a lab can see which analyses are affected.
    files = dir(fullfile(folders{k},'*.m'));
    names = {}; callers = {};
    for f = 1:numel(files)
        code = fileread(fullfile(files(f).folder,files(f).name));
        code = regexprep(code,'(?m)^\s*%[^\r\n]*','');
        calls = regexp(code, pattern, 'match');
        if isempty(calls), continue; end
        local = regexp(code,'(?m)^\s*function\s+(?:[^=\r\n]*=\s*)?(\w+)\s*\(', 'tokens');
        local = cellfun(@(m) m{1}, local, 'UniformOutput', false);
        absent = setdiff(calls, local);
        if ~isempty(absent)
            names = union(names, absent);
            callers{end+1} = files(f).name; %#ok<AGROW>
        end
    end
    if ~isempty(names)
        [~,component] = fileparts(folders{k});
        rows(end+1,:) = {['NF helper references in ' component], ...
            nf_probe_result(names), ['Referenced by ' strjoin(callers, ', ')]}; %#ok<AGROW>
    end
end
end
