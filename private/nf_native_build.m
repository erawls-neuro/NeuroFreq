function [native,rows] = nf_native_build(root,mode)
% Build only NF's named entry points in the two source directories below.
% Vendor packages manage their own binaries. An unsuccessful target leaves
% its MATLAB implementation available and never promotes a partial build.
native = struct('abi',nf_native_abi(),'mode',mode,'targets',struct([]));
rows = cell(0,3);
folders = {fullfile('phaseamplitudecoupling','hnf_phaseamplitudecoupling'), ...
    fullfile('timefrequency','hnf_timefrequency')};
toolchains = struct();
cache = fullfile(root,'dependencies','_nf_local','native');
for i = 1:numel(folders)
    listing = dir(fullfile(root,folders{i},'hnf_*_mex.*'));
    for j = 1:numel(listing)
        filename = listing(j).name;
        if listing(j).isdir || isempty(regexp(filename,'^hnf_[A-Za-z0-9_]+_mex\.(c|cpp|cu)$','once'))
            continue
        end
        [~,name,ext] = fileparts(filename);
        status = 'skipped'; detail = 'MATLAB implementation selected.';
        if strcmp(mode,'auto')
            try
                family = ext(2:end);
                if ~isfield(toolchains,family)
                    toolchains.(family) = compiler_info(family);
                end
                tool = toolchains.(family);
                if ~isempty(tool.problem), error('NeuroFreq:NativeCompiler','%s',tool.problem); end
                [record,status] = install_one(root,folders{i},filename,cache,native.abi,tool.signature);
                if ~isempty(native.targets) && any(strcmp({native.targets.name},name))
                    error('NeuroFreq:NativeName','Duplicate MEX entry point %s.',name);
                end
                native.targets(end+1) = record;
                detail = record.binary;
            catch err
                status = 'fallback'; detail = err.message;
            end
        end
        fprintf('  %s: %s. %s\n',name,status,detail);
        rows(end+1,:) = {name,status,detail}; %#ok<AGROW>
    end
end
end

function info = compiler_info(family)
% Record selected compiler and actual optimization/ABI options. Reusing a
% binary for the same interpreter does not require checking out a license.
info = struct('signature','','problem','');
try
    if exist('OCTAVE_VERSION','builtin')
        if strcmp(family,'cu')
            error('Octave CUDA builds are unavailable.');
        end
        vars = {'CC','CXX','CPPFLAGS','CFLAGS','CXXFLAGS','LDFLAGS','DL_LDFLAGS'};
        values = cell(size(vars));
        for i=1:numel(vars)
            [values{i},status] = mkoctfile('-p',vars{i});
            if status~=0, error('Could not identify the Octave native toolchain.'); end
        end
        % The compiler used to build Octave and its development headers is
        % part of the ABI. Commands/flags can be overridden by its environment.
        compilerIndex = 1 + strcmp(family,'cpp');
        compiler = values{compilerIndex};
        % Ask a simple compiler executable directly (no shell evaluation).
        % Also keep the full mkoctfile configuration for wrapped compilers.
        info.signature = strjoin([{'mkoctfile --mex -O',compiler_banner(compiler)},values],'|');
    else
        language = 'C';
        if ~strcmp(family,'c'), language='C++'; end
        compiler = mex.getCompilerConfigurations(language,'Selected');
        if isempty(compiler)
            error('No selected %s compiler. Configure one with mex -setup %s.',language,language);
        end
        signature = {compiler.Name, compiler.Version, compiler.Location, compiler.LinkerName, ...
            compiler.LinkerVersion, '-O -R2018a'};
        if isfile(compiler.MexOpt)
            signature{end+1} = nf_native_digest({compiler.MexOpt},'mex options');
        end
        % Version and Location may be empty on Unix. Capture the configured
        % command, its flags, and the compiler's own version where available.
        details = compiler.Details;
        signature{end+1} = evalc('disp(details)');
        if isunix && isprop(details,'CompilerExecutable')
            signature{end+1} = compiler_banner(details.CompilerExecutable);
        end
        vars = {'CC','CXX','CFLAGS','CXXFLAGS','LDFLAGS','COMPFLAGS','LINKFLAGS', ...
            'NVCCFLAGS','MW_NVCC_PATH','CUDA_PATH'};
        signature = [signature,cellfun(@getenv,vars,'UniformOutput',false)];
        if strcmp(family,'cu')
            if isempty(which('mexcuda'))
                error('mexcuda is unavailable for this MATLAB installation.');
            end
            % MATLAB ships nvcc with the release recorded in the ABI. When
            % a separate toolkit is selected, include its compiler version.
            override = getenv('MW_NVCC_PATH');
            if ~isempty(override)
                if isfolder(override), override = fullfile(override,'nvcc'); end
                signature{end+1} = compiler_banner(override);
            end
        end
        info.signature = strjoin(cellfun(@char,signature,'UniformOutput',false),'|');
    end
catch err
    info.problem = err.message;
end
end

function banner = compiler_banner(executable)
% Ask a simple compiler executable directly. Wrapped shell commands and
% placeholders remain part of the key as text and are never evaluated.
banner = char(executable);
if isempty(regexp(banner,'^[A-Za-z0-9_./:+-]+$','once')) && ~isfile(banner)
    return
end
args = javaArray('java.lang.String',2);
args(1) = javaObject('java.lang.String',banner);
args(2) = javaObject('java.lang.String','--version');
process = javaObject('java.lang.ProcessBuilder',args);
child = process.redirectErrorStream(true).start();
reader = javaObject('java.io.BufferedReader',javaObject('java.io.InputStreamReader',child.getInputStream()));
lines = {};
line = reader.readLine();
while ~isempty(line)
    lines{end+1} = char(line); %#ok<AGROW>
    line = reader.readLine();
end
reader.close();
child.waitFor();
banner = strjoin([{banner},lines],' ');
end

function [item,status] = install_one(root,folder,filename,cache,abi,tool)
source = fullfile(root,folder,filename);
inputs = nf_native_inputs(source);
sourceHash = nf_native_digest(inputs,'NF source v1');
key = nf_native_digest({},strjoin({abi,tool,sourceHash},'|'));
[~,name,ext] = fileparts(source);
bin = [name '.' mexext];
item = struct('name',name,'source',fullfile(folder,filename),'sourceHash',sourceHash, ...
    'binary','','binaryHash','');
if ~isfolder(cache)
    [ok,msg] = mkdir(cache);
    if ~ok, error('NeuroFreq:NativeCache','%s',msg); end
end
if ~strcmp(nf_canonical(cache),cache)
    error('NeuroFreq:NativeCache','Native cache resolves outside its intended location.');
end
found = dir(fullfile(cache,[key '-*']));
for i=1:numel(found)
    if ~found(i).isdir, continue; end
    location = fullfile(cache,found(i).name);
    binary = fullfile(location,bin);
    receipt = fullfile(location,'sha256.txt');
    directory = javaObject('java.io.File',location);
    members = cellstr(directory.list());
    if ~isequal(sort(members(:)),sort({bin;'sha256.txt'})), continue; end
    if ~isfile(binary) || ~isfile(receipt) || ~strcmp(nf_canonical(binary),binary), continue; end
    try
        actual = nf_native_digest({binary},'NF binary v1');
        valid = strcmp(strtrim(fileread(receipt)),actual);
    catch
        valid = false;
    end
    if valid
        item.binary = binary(numel(root)+2:end);
        item.binaryHash = actual;
        status = 'cached';
        return
    end
end

stage = tempname(cache);
[ok,msg] = mkdir(stage);
if ~ok, error('NeuroFreq:NativeStage','%s',msg); end
cleanup = onCleanup(@() remove_stage(stage)); %#ok<NASGU>
% Compile from an isolated working directory, so object files and temporary
% linker output never land in the user's working directory or the NF sources.
for i=1:numel(inputs), copyfile(inputs{i},stage); end
if ~strcmp(nf_native_digest(nf_native_inputs(fullfile(stage,filename)),'NF source v1'),sourceHash)
    error('NeuroFreq:NativeChanged','Source changed while staging. Run nf_install again.');
end
compile_in(stage,filename,ext);
binary = fullfile(stage,bin);
if ~isfile(binary), error('NeuroFreq:NativeBuild','Compiler did not produce %s.',bin); end
item.binaryHash = nf_native_digest({binary},'NF binary v1');
% Publish only the binary and its checksum. Keep immutable generations so a
% failed or repeated install never overwrites a binary loaded by MATLAB.
publish = fullfile(stage,'out'); mkdir(publish);
movefile(binary,fullfile(publish,bin));
fid = fopen(fullfile(publish,'sha256.txt'),'wt');
if fid<0, error('NeuroFreq:NativeBuild','Could not write build receipt.'); end
fprintf(fid,'%s\n',item.binaryHash); fclose(fid);
[~,uniqueName] = fileparts(stage);
final = fullfile(cache,[key '-' uniqueName]);
from = javaObject('java.io.File',publish);
to = javaObject('java.io.File',final);
options = javaArray('java.nio.file.CopyOption',1);
options(1) = javaMethod('valueOf','java.nio.file.StandardCopyOption','ATOMIC_MOVE');
javaMethod('move','java.nio.file.Files',from.toPath(),to.toPath(),options);
binary = fullfile(final,bin);
item.binary = binary(numel(root)+2:end);
status = 'built';
end

function compile_in(stage,filename,ext)
original=pwd; restore=onCleanup(@() cd(original)); %#ok<NASGU>
cd(stage);
if exist('OCTAVE_VERSION','builtin')
    [message,status] = mkoctfile('--mex','-O',filename);
    if status~=0, error('NeuroFreq:NativeBuild','mkoctfile failed. %s',message); end
elseif strcmp(ext,'.cu')
    mexcuda('-O','-R2018a',filename);
else
    mex('-O','-R2018a',filename);
end
end

function remove_stage(folder)
if isfolder(folder), rmdir(folder,'s'); end
end
