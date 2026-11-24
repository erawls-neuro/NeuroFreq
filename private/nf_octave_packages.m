function rows = nf_octave_packages()
% Use installed Octave packages within the caller's path scope. No download
% or permanent package activation. pkg also loads declared dependencies.
installed = pkg('list');
names = cellfun(@(p) p.name,installed,'UniformOutput',false);
% These packages define competing table and categorical classes. Load only
% one provider. Current Octave statistics depends on the datatypes package.
if any(strcmp(names,'datatypes')), tables = 'datatypes'; else, tables = 'tablicious'; end
packages = {'signal','statistics','optim','io',tables};
uses = {'Filter design, analytic signals and spectra', ...
    'Probability distributions and multivariate statistics', ...
    'Constrained optimization', 'Spreadsheet import', 'Table and categorical data'};
rows = cell(numel(packages),3);
% Statistics deliberately shadows some core reductions. Keep the user's
% warning preference outside package registration, including after errors.
shadow = warning('query','Octave:shadowed-function');
restore = onCleanup(@() warning(shadow.state,shadow.identifier)); %#ok<NASGU>
warning('off','Octave:shadowed-function');
for k = 1:numel(packages)
    name = packages{k};
    index = find(strcmp(name,names),1);
    status = 'Not installed.';
    if ~isempty(index)
        try
            if ~installed{index}.loaded, pkg('load',name); end
            status = ['Available (' installed{index}.version ').'];
        catch problem
            status = ['Load failed: ' problem.message];
            warning('NeuroFreq:OctavePackage','%s: %s',name,status);
        end
    end
    rows(k,:) = {['Octave package ' name],status,uses{k}};
end
end
