function [tab, behMap] = nf_unpackdata(TF, varargin)
% NF_UNPACKDATA    Flatten a NeuroFreq TF structure into a MATLAB table
%
% Syntax
%   [TAB, BEHMAP] = NF_UNPACKDATA(TF)                   % returns wide table
%   [TAB, BEHMAP] = NF_UNPACKDATA(TF,'Mode','long',..)  % returns long table
%
% Mode
%   'wide' (default)  row = {time, trial}
%                     column = {measure, channel, frequency}
%   'long'            row = {channel, frequency, time, trial}
%                     plus one column per selected measure
%
% Measure flags (logical; defaults in brackets)
%   'IncludePower'   TF.power                [true]
%   'IncludePhase'   TF.phase                [false]
%   'IncludeOsc'     TF.SPRiNT.osc_power     [false]
%   'IncludeAper'    TF.SPRiNT.ap_power      [false]
%   'IncludeERP'     TF.erprem.erppow        [false]
%   'IncludeERPRem'  TF.erprem.erprempow     [false]
%
% Index slices ([] keeps all; numeric or logical)
%   'ChanIdx' , 'FreqIdx' , 'TimeIdx' , 'TrialIdx'
%
% behMap
%   A struct mapping each string behavioral variable in TF.behavior
%   to its categorical codes (empty if no string behaviors present).
%
% Rules
%   * All tensors must be size [nChan nFreq nTime nTrial].
%   * Missing trial dimension (resting-state) is allowed and treated as 1.
%   * Behavioral variables must be scalar per trial.
%
% Errors
%   nf_unpackdata:DimMismatch   size mismatch in a requested tensor
%   nf_unpackdata:BehaviorShape behavioral variable not scalar per trial
%   nf_unpackdata:MissingField  expected field absent
%   nf_unpackdata:*Bound        illegal index
%   nf_unpackdata:*Len          logical index wrong length
%


% option parsing -----------------------------------------------------------
tfFlag = @(b) islogical(b) && isscalar(b);
tfIdx  = @(v) isnumeric(v) || islogical(v);

ip = inputParser;
ip.addParameter('Mode','wide',@(s)any(strcmpi(s,{'wide','long'})));
ip.addParameter('IncludePower',  true ,tfFlag);
ip.addParameter('IncludePhase',  false,tfFlag);
ip.addParameter('IncludeOsc',    false,tfFlag);
ip.addParameter('IncludeAper',   false,tfFlag);
ip.addParameter('IncludeERP',    false,tfFlag);
ip.addParameter('IncludeERPRem', false,tfFlag);
ip.addParameter('ChanIdx',  [],tfIdx);
ip.addParameter('FreqIdx',  [],tfIdx);
ip.addParameter('TimeIdx',  [],tfIdx);
ip.addParameter('TrialIdx', [],tfIdx);
ip.parse(varargin{:});
o        = ip.Results;
wideMode = strcmpi(o.Mode,'wide');

% dimensions ---------------------------------------------------------------
nC = TF.nsensor;                                    % channels (mandatory)
nF = numel(TF.freqs);                               % frequencies
nT = numel(TF.times);                               % time points

% trials may be absent for resting-state; fall back to tensor size
if isfield(TF,'ntrls') && ~isempty(TF.ntrls)
    nR = TF.ntrls;
else
    nR = size(getFirstTensor(TF),4);                % at least 1
    if nR == 0, nR = 1; end
end

iC = makeIndex(o.ChanIdx , nC , 'ChanIdx');
iF = makeIndex(o.FreqIdx , nF , 'FreqIdx');
iT = makeIndex(o.TimeIdx , nT , 'TimeIdx');
iR = makeIndex(o.TrialIdx, nR , 'TrialIdx');

labC = string({TF.chanlocs(iC).labels});
valF = TF.freqs(iF);

% coordinate columns -------------------------------------------------------
if wideMode
    tCol = repmat(TF.times(iT).', numel(iR), 1);
    rMat = kron(iR.', ones(1,numel(iT))).';
    coordTime  = tCol(:);
    coordTrial = rMat(:);
    if max(coordTrial) <= intmax('uint16')
        coordTrial = uint16(coordTrial);
    else
        coordTrial = uint32(coordTrial);
    end
    tab = table(coordTime, coordTrial, 'VariableNames',{'time','trial'});
else
    nc = numel(iC); nf = numel(iF); nt = numel(iT); nr = numel(iR);
    chanCol = kron(ones(nf*nt*nr,1), labC(:));
    frBase  = kron(ones(nt*nr,1), valF.');
    freqCol = repmat(frBase, nc, 1);
    tBase   = repmat(TF.times(iT).', nr, 1);
    timeCol = repmat(tBase(:), nc*nf, 1);
    rBase   = repelem(iR(:).', nt, 1);
    trialCol= repmat(rBase(:), nc*nf, 1);
    if max(trialCol) <= intmax('uint16')
        trialCol = uint16(trialCol);
    else
        trialCol = uint32(trialCol);
    end
    tab = table(categorical(chanCol), freqCol, timeCol, trialCol,...
                'VariableNames',{'channel','frequency','time','trial'});
end

% selection map and helpers -----------------------------------------------
need = struct('power',o.IncludePower,'phase',o.IncludePhase,...
              'osc',o.IncludeOsc,'aper',o.IncludeAper,...
              'erp',o.IncludeERP,'erpr',o.IncludeERPRem);

slice4  = @(X) X(iC,iF,iT,iR);
vec     = @(X) reshape(slice4(X),[],1);
addWide = @(X,pfx) makeWideCols(slice4(X),pfx,labC,valF,tab);

% numeric blocks -----------------------------------------------------------
if need.power && isfield(TF,'power')
    chkSize(TF.power,'power',nC,nF,nT,nR);
    if wideMode, tab = [tab addWide(TF.power,'Pwr')]; else, tab.power = vec(TF.power); end
end
if need.phase && isfield(TF,'phase')
    chkSize(TF.phase,'phase',nC,nF,nT,nR);
    if wideMode, tab = [tab addWide(TF.phase,'Phs')]; else, tab.phase = vec(TF.phase); end
end
if need.osc
    chkFld(TF,'SPRiNT'); chkFld(TF.SPRiNT,'osc_power');
    chkSize(TF.SPRiNT.osc_power,'osc_power',nC,nF,nT,nR);
    if wideMode, tab = [tab addWide(TF.SPRiNT.osc_power,'Osc')]; else, tab.osc_power = vec(TF.SPRiNT.osc_power); end
end
if need.aper
    chkFld(TF,'SPRiNT'); chkFld(TF.SPRiNT,'ap_power');
    chkSize(TF.SPRiNT.ap_power,'ap_power',nC,nF,nT,nR);
    if wideMode, tab = [tab addWide(TF.SPRiNT.ap_power,'Aper')]; else, tab.ap_power = vec(TF.SPRiNT.ap_power); end
end
if need.erp
    chkFld(TF,'erprem'); chkFld(TF.erprem,'erppow');
    chkSize(TF.erprem.erppow,'erppow',nC,nF,nT,nR);
    if wideMode, tab = [tab addWide(TF.erprem.erppow,'ERP')]; else, tab.erp_power = vec(TF.erprem.erppow); end
end
if need.erpr
    chkFld(TF,'erprem'); chkFld(TF.erprem,'erprempow');
    chkSize(TF.erprem.erprempow,'erprempow',nC,nF,nT,nR);
    if wideMode, tab = [tab addWide(TF.erprem.erprempow,'ERPr')]; else, tab.erprem_power = vec(TF.erprem.erprempow); end
end

% behavioral variables -----------------------------------------------------
behMap = struct();
if isfield(TF,'behavior') && ~isempty(TF.behavior)
    fn = fieldnames(TF.behavior);
    for k = 1:numel(fn)
        nm = fn{k};
        vals = [TF.behavior.(nm)];
        if ~isscalar(vals(1))
            error('nf_unpackdata:BehaviorShape','Behavior %s not scalar per trial.',nm);
        end
        vals = vals(iR);
        if ischar(vals) || isstring(vals) || iscellstr(vals)
            cat = categorical(vals);
            behMap.(nm).categories = categories(cat);
            vals = double(cat);
        end
        if wideMode
            tab.(nm) = repelem(vals(:), numel(iT));
        else
            rep = numel(iC)*numel(iF);
            tab.(nm) = kron(vals(:), ones(numel(iT)*rep,1));
        end
    end
end
end

% helper: first available tensor to infer trial dimension
function T = getFirstTensor(TF)
fields = {'power','phase'};
for f = 1:numel(fields)
    if isfield(TF,fields{f}); T = TF.(fields{f}); return; end
end
error('nf_unpackdata:NoTensor','No numeric tensors found to infer size.');
end

% helper: index builder
function idx = makeIndex(u,dim,tag)
if isempty(u), idx = 1:dim; return, end
if islogical(u)
    if numel(u)~=dim, error('nf_unpackdata:%sLen',tag); end
    idx = find(u);
else
    if any(u<1|u>dim), error('nf_unpackdata:%sBound',tag); end
    idx = u(:).';
end
end

% helper: wide column construction
function T = makeWideCols(B,prefix,chanLabs,freqVals,baseTab)
[mC,mF,~,~] = size(B);
M = reshape(B,mC*mF,[]).';
cChan = repmat(chanLabs,1,mF).';
cFreq = repelem(freqVals,mC).';
raw   = compose('%s_%s_%.4gHz',prefix,cChan,cFreq);
names = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(raw), baseTab.Properties.VariableNames);
T = array2table(M,'VariableNames',names);
end

% helper: tensor size check
function chkSize(A,tag,a,b,c,d)
sz = size(A); sz(end+1:4)=1; % pad missing dims (trials may be absent)
if ~isequal(sz,[a b c d])
    error('nf_unpackdata:DimMismatch','Matrix %s wrong size.',tag);
end
end

% helper: field presence
function chkFld(S,f)
if ~isfield(S,f), error('nf_unpackdata:MissingField','Missing field %s.',f); end
end
