function report = nf_configure(root, report)
% Configure the copy that will actually run. Keep machine-specific state in
% one untracked folder. Third-party packages and MATLAB paths never belong
% in the reference inventory used to verify a GitHub release.
cache = fullfile(root, 'dependencies', '_nf_local');
if ~strcmp(nf_canonical(cache), cache)
    error('NeuroFreq:LinkedCache', ...
        'The installer cache resolves to another location. Inspect %s before installing.', cache);
end
if ~isfolder(cache)
    [ok, message] = mkdir(cache);
    if ~ok, error('NeuroFreq:Cache', 'Cannot create %s. %s', cache, message); end
end

fprintf('\nEEGLAB\n');
[eeglabRoot, privateEEGLAB, report.eeglab] = nf_eeglab(root);
fprintf('%s\n', report.eeglab);
fprintf('\nNeuroFreq dependencies and EEGLAB plugins\n');
[report.dependencies, paths] = nf_dependencies(root, eeglabRoot, privateEEGLAB);
fprintf('\nMATLAB feature availability\n');
report.features = nf_check_features(paths, root, eeglabRoot);

% nf_run temporarily adds the chosen providers before dispatching a call.
% Store local paths relative to this checkout so a moved NF folder still works.
% Plugin and dependency directories precede EEGLAB's general function folders.
eegPaths = {};
if ~isempty(eeglabRoot), eegPaths = nf_eeglab_paths(eeglabRoot); end
paths = unique([paths(:); eegPaths(:)], 'stable');
config.scopePaths = cellfun(@(p) nf_portable_path(p, root), paths, ...
    'UniformOutput', false);
config.eeglabRoot = nf_portable_path(eeglabRoot, root);
config.privateEEGLAB = privateEEGLAB;

% A failed write must leave the last usable configuration intact. A temporary
% file in the same folder can be renamed over runtime.mat in one operation.
staged = [tempname(cache) '.mat'];
cleanup = onCleanup(@() nf_remove_temp_file(staged)); %#ok<NASGU>
save(staged, 'config', '-v7');
source = javaObject('java.io.File', staged);
target = javaObject('java.io.File', fullfile(cache, 'runtime.mat'));
options = javaArray('java.nio.file.CopyOption', 2);
options(1) = javaMethod('valueOf', 'java.nio.file.StandardCopyOption', 'REPLACE_EXISTING');
options(2) = javaMethod('valueOf', 'java.nio.file.StandardCopyOption', 'ATOMIC_MOVE');
javaMethod('move', 'java.nio.file.Files', source.toPath(), target.toPath(), options);

report.configured = true;
fprintf('\nNeuroFreq configuration saved. Dependency results:\n');
for i = 1:size(report.dependencies, 1)
    fprintf('  %s: %s\n', report.dependencies{i, 1}, report.dependencies{i, 2});
end
fprintf(['\nRun an NF function with its configured dependencies:\n' ...
    '  EEG = nf_run(''nf_filter'', EEG);\n' ...
    'Run a pipeline containing ordinary nf_* calls:\n' ...
    '  nf_run(@my_analysis);\n']);
end

function p = nf_portable_path(p, root)
if isempty(p), p = ''; return; end
p = nf_canonical(p);
if nf_under(p, root)
    if numel(p) == numel(root)
        p = '.';
    else
        p = p(numel(root)+2:end);
    end
end
end

function nf_remove_temp_file(file)
if isfile(file), delete(file); end
end
