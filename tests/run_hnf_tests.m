function report = run_hnf_tests
% RUN_HNF_TESTS  Run all supplied base-MATLAB numerical contract tests.
%   R=run_hnf_tests; assert(all([R.passed]));
%   Captured console/warning text fails a test, preserving the quiet API.
%   Restores the caller's MATLAB path and random stream on exit.
folder = fileparts(mfilename('fullpath'));
root = fileparts(folder); oldPath = path; oldRng = rng;
cleanup = onCleanup(@() restore_environment(oldPath,oldRng)); %#ok<NASGU>
addpath(folder,'-begin'); addpath(fullfile(root,'helpers','robustarithmetic'),'-begin');
files = dir(fullfile(folder,'test_hnf_*.m'));
[~,order] = sort({files.name}); files=files(order);
report = repmat(struct('name','','passed',false,'seconds',NaN, ...
    'errorIdentifier','','message','','matlabVersion',version),numel(files),1);
for k=1:numel(files)
    [~,testName] = fileparts(files(k).name); report(k).name = testName;
    started = tic;
    try
        consoleText = evalc('feval(testName);');
        assert(isempty(strtrim(consoleText)),'hnf:test:ConsoleOutput', ...
            'A test or low-level function produced console/warning output: %s',consoleText);
        report(k).passed = true;
    catch ME
        report(k).errorIdentifier = ME.identifier;
        report(k).message = ME.message;
    end
    report(k).seconds = toc(started);
end
end
function restore_environment(oldPath,oldRng)
path(oldPath); rng(oldRng);
end
