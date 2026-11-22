function test_nf_install()
% TEST_NF_INSTALL  Offline regression tests for update and path isolation.
% From the repository root, run: addpath(fullfile(pwd,'tests')); test_nf_install
% All changes take place in a temporary directory. The private fixture
% functions below replace GitHub and the desktop dialogs only in that copy.
repo = fileparts(fileparts(mfilename('fullpath')));
previousPath = path;
previousFolder = pwd;
work = char(javaObject('java.io.File',tempname).getCanonicalPath());
mkdir(work);
cleanup = onCleanup(@() finish(previousFolder, previousPath, work)); %#ok<NASGU>
checkUpdate(repo, work);
checkRuntime(repo, work);
fprintf('NeuroFreq installer regression tests passed.\n');
end

function checkUpdate(repo, work)
root = fullfile(work, 'NeuroFreq');
private = fullfile(root, 'private');
mkdir(private);
copyfile(fullfile(repo, 'nf_install.m'), root);
writeText(fullfile(root, 'personal.m'), 'my local edit');
% Read fixture state from the common parent, even after a directory rename.
stateFile = 'fullfile(fileparts(fileparts(fileparts(mfilename(''fullpath'')))), ''state.mat'')';
writeText(fullfile(private,'webread.m'), {
    'function value = webread(url,varargin)'
    ['state = load(' stateFile ');']
    'if ~isempty(strfind(url,''/commits'')), value=state.head; else, value=state.tree; end'
    'end'});
writeText(fullfile(private,'websave.m'), {
    'function dest = websave(dest,varargin)'
    ['state = load(' stateFile '); copyfile(state.archive,dest);']
    'end'});
writeText(fullfile(private,'questdlg.m'), {
    'function value = questdlg(varargin)'
    ['state = load(' stateFile '); value=state.answer;']
    'end'});
writeText(fullfile(private,'usejava.m'), {'function yes=usejava(varargin)', 'yes=true;', 'end'});
writeText(fullfile(private,'weboptions.m'), {'function opt=weboptions(varargin)', 'opt=struct();', 'end'});
writeText(fullfile(private,'tempdir.m'), {'function value=tempdir', ...
    'value=fileparts(fileparts(fileparts(mfilename(''fullpath''))));', 'end'});

% Known Git blob hash for the two lines "function nf_install\nend\n".
entry = struct('type','blob','mode','100644','path','nf_install.m', ...
    'sha','dbfa66c7e887c84ca5167bc1fa3f8297bc3edb34');
tree = struct('sha',repmat('a',1,40),'truncated',false,'tree',entry);
head = struct('sha',repmat('b',1,40),'commit',struct('tree',struct('sha',tree.sha)));
archive = fullfile(work, 'release.zip');
answer = 'Keep this copy';
save(fullfile(work,'state.mat'), 'tree', 'head', 'archive', 'answer');

cd(root);
clear nf_install
report = nf_install();
assert(strcmp(report.repository,'kept') && isfile(fullfile(root,'personal.m')));
assert(~isfolder([root '_deprecated']));

% A corrupt archive must leave the original in place after approval.
answer = 'Update';
save(fullfile(work,'state.mat'), 'tree', 'head', 'archive', 'answer');
release = fullfile(work,'release');
mkdir(release);
writeText(fullfile(release,'nf_install.m'), 'corrupt');
zip(archive, {'release'}, work);
report = nf_install();
assert(strcmp(report.repository,'unverified') && ~isfolder([root '_deprecated']));
assert(strcmp(fileread(fullfile(root,'personal.m')), 'my local edit'));

% A valid archive keeps both the user's modified tree and an older backup.
delete(archive);
writeText(fullfile(release,'nf_install.m'), {'function nf_install','end'});
zip(archive, {'release'}, work);
mkdir([root '_deprecated']);
writeText(fullfile([root '_deprecated'],'old.txt'), 'previous backup');
report = nf_install();
assert(strcmp(report.repository,'updated'));
assert(strcmp(fileread(fullfile([root '_deprecated'],'old.txt')), 'previous backup'));
assert(strcmp(fileread(fullfile([root '_deprecated_2'],'personal.m')), 'my local edit'));
assert(~isfile(fullfile(root,'personal.m')));
assert(strcmp(fileread(fullfile(root,'nf_install.m')), sprintf('function nf_install\nend\n')));
cd(work);
end

function checkRuntime(repo, work)
root = fullfile(work,'runtime');
mkdir(root);
copyfile(fullfile(repo,'nf_run.m'),root);
copyfile(fullfile(repo,'private'),fullfile(root,'private'));
cache = fullfile(root,'dependencies','_nf_local');
local = fullfile(cache,'EEGLAB'); external = fullfile(work,'external');
mkdir(local); mkdir(external);
writeText(fullfile(local,'eeglab.m'), {'function value=eeglab','value=1;','end'});
writeText(fullfile(external,'eeglab.m'), {'function value=eeglab','value=2;','end'});
mkdir(fullfile(root,'preprocessing'));
writeText(fullfile(root,'preprocessing','nf_fixture.m'), {
    'function [value,location]=nf_fixture(fail)'
    'value=eeglab(); location=which(''eeglab'');'
    'if fail, rmpath(fileparts(location)); error(''NFTest:Expected'',''test error''); end'
    'end'});
config = struct('scopePaths',{{}}, 'eeglabRoot',fullfile('dependencies','_nf_local','EEGLAB'), ...
    'privateEEGLAB',true);
save(fullfile(cache,'runtime.mat'),'config');
addpath(root,external,'-begin');
clear nf_run
original = path;
assert(eeglab()==2);
[value,location] = nf_run('nf_fixture',false);
assert(value==1 && strcmp(location,fullfile(local,'eeglab.m')) && strcmp(original,path));
try
    nf_run('nf_fixture',true);
    error('NFTest:MissingError','Expected the fixture to throw.');
catch failure
    assert(strcmp(failure.identifier,'NFTest:Expected'),failure.message);
end
assert(strcmp(original,path) && eeglab()==2);
end

function writeText(file, text)
% Write binary so fixture hashes are identical on Windows, macOS and Linux.
fid = fopen(file,'wb');
assert(fid >= 0);
closeFile = onCleanup(@() fclose(fid)); %#ok<NASGU>
if iscell(text), text = sprintf('%s\n',text{:}); end
fwrite(fid,text,'char');
end

function finish(previousFolder, previousPath, work)
cd(previousFolder);
path(previousPath);
clear nf_install nf_run
rmdir(work,'s');
end
