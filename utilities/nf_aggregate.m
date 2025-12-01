function tfagg = nf_aggregate( )
% NF_AGGREGATE    Aggregates averaged TF datasets into a multi-subject set.
%
%
% OUTPUT
% ------
% TF - multi-subject TF structure
%
% INPUT
% -----
% no input - user select
%
%


%get files
[nm,pt] = uigetfile('*.*',  'All Files (*.*)','MultiSelect','on');

%get all data
for i=1:numel(nm)
    %get tf structure
    fi=load([pt filesep nm{i}],'-mat');
    f = fieldnames(fi);
    tf = fi.(f{1});
    %set data sizes
    if i==1
        tfagg.power = zeros( numel(nm),length({tf.chanlocs.labels}),numel(tf.freqs),numel(tf.times),tf.conds );
        if isfield(tf,'phase')
            tfagg.phase = zeros( numel(nm),length({tf.chanlocs.labels}),numel(tf.freqs),numel(tf.times),tf.conds );
        end
        tfagg.nsubjects = numel(nm);
        tfagg.freqs = tf.freqs;
        tfagg.times = tf.times;
        tfagg.Fs = tf.Fs;
        tfagg.conds = tf.conds;
        if tf.ntrls ~= 1
            error('must be averaged TF');
        end
        tfagg.scale = tf.scale;
        tfagg.chanlocs = tf.chanlocs;
        if isfield(tf,'SPRiNT')
            tfagg.SPRiNT.ap_power = zeros( numel(nm), numel({tfagg.chanlocs.labels}), numel(tfagg.freqs), numel(tfagg.times), tf.conds );
            tfagg.SPRiNT.osc_power = zeros( numel(nm), numel({tfagg.chanlocs.labels}), numel(tfagg.freqs), numel(tfagg.times), tf.conds );
            tfagg.SPRiNT.rsquare = zeros( numel(nm), numel({tfagg.chanlocs.labels}), numel(tfagg.times), tf.conds );
        end
        if isfield(tf,'cpm')
            tfagg.cpm.PLmean = zeros( numel(nm), numel({tfagg.chanlocs.labels}), numel(tfagg.freqs), numel(tfagg.times), tf.conds );
            tfagg.cpm.NPLmean = zeros( numel(nm), numel({tfagg.chanlocs.labels}), numel(tfagg.freqs), numel(tfagg.times), tf.conds );
            tfagg.cpm.PLvar = zeros( numel(nm), numel({tfagg.chanlocs.labels}), numel(tfagg.freqs), numel(tfagg.times), tf.conds );
            tfagg.cpm.NPLvar = zeros( numel(nm), numel({tfagg.chanlocs.labels}), numel(tfagg.freqs), numel(tfagg.times), tf.conds );
            tfagg.cpm.alpha = zeros( numel(nm), numel({tfagg.chanlocs.labels}), numel(tfagg.freqs), numel(tfagg.times), tf.conds );
        end
        if isfield(tf,'erprem')
            tfagg.erprem.erprempow = zeros( numel(nm), numel({tfagg.chanlocs.labels}), numel(tfagg.freqs), numel(tfagg.times), tf.conds );
            tfagg.erprem.erppow = zeros( numel(nm), numel({tfagg.chanlocs.labels}), numel(tfagg.freqs), numel(tfagg.times), tf.conds );
            tfagg.erprem.erpremphase = zeros( numel(nm), numel({tfagg.chanlocs.labels}), numel(tfagg.freqs), numel(tfagg.times), tf.conds );
        end
        if isfield(tf,'behavior')
            if length( tf.behavior ) ~= tf.conds
                error('behavior and data are not equivalent trial lengths?');
            end
            tfagg.behavior = tf.behavior;
        end
        subvec = [];
        condvec = [];
    else
        if ~all(tf.freqs==tfagg.freqs)==1 || ~all(tf.times==tfagg.times)==1 || tf.Fs~=tfagg.Fs || tf.ntrls~=1 || strcmp(tf.scale,tfagg.scale)~=1 || ~all([tfagg.chanlocs.labels]==[tf.chanlocs.labels])==1
            error('dimension mismatch!');
        end
        if isfield(tfagg,'behavior')
            if length( tf.behavior )~=tfagg.conds
                error('dimension mismatch!');
            end
        end
    end
    tfagg.power(i,:,:,:,:) = tf.power;
    if isfield(tfagg,'phase')
        tfagg.phase(i,:,:,:,:) = tf.phase;
    end
    if isfield(tfagg,'SPRiNT')
        tfagg.SPRiNT.ap_power(i,:,:,:,:) = tf.SPRiNT.ap_power;
        tfagg.SPRiNT.osc_power(i,:,:,:,:) = tf.SPRiNT.osc_power;
        tfagg.SPRiNT.rsquare(i,:,:,:) = tf.SPRiNT.rsquare;
    end
    if isfield(tf,'cpm')
        tfagg.cpm.PLmean(i,:,:,:,:) = tf.cpm.PLmean;
        tfagg.cpm.NPLmean(i,:,:,:,:) = tf.cpm.NPLmean;
        tfagg.cpm.PLvar(i,:,:,:,:) = tf.cpm.PLvar;
        tfagg.cpm.NPLvar(I,:,:,:,:) = tf.cpm.NPLvar;
        tfagg.cpm.alpha(i,:,:,:,:) = tf.cpm.alpha;
    end
    if isfield(tf,'erprem')
        tfagg.erprem.erprempow(i,:,:,:,:) = tf.erprem.erprempow;
        tfagg.erprem.erppow(i,:,:,:,:) = tf.erprem.erppow;
        tfagg.erprem.erpremphase(i,:,:,:,:) = tf.erprem.erpremphase;
    end
    if isfield(tf,'behavior') && i>1
        tfagg.behavior = [tfagg.behavior tf.behavior];
    end
    subvec = [subvec repmat(i,[1 tfagg.conds])];
    condvec = [condvec 1:tfagg.conds];
end
for i=1:length(subvec)
    tfagg.behavior(i).subID = subvec(i);
    tfagg.behavior(i).condition = condvec(i);
end  
tfagg.subid = nm;






