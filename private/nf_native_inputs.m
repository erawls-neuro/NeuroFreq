function files = nf_native_inputs(source)
% Native kernels are stand-alone translation units. A shared header belongs
% beside its kernels. Include all sibling headers in the build fingerprint.
folder = fileparts(source);
headers = [dir(fullfile(folder,'*.h')); dir(fullfile(folder,'*.hpp')); dir(fullfile(folder,'*.cuh'))];
if isempty(headers), files = {source}; return; end
names = sort({headers(~[headers.isdir]).name});
files = [{source}, cellfun(@(n) fullfile(folder,n),names,'UniformOutput',false)];
end
