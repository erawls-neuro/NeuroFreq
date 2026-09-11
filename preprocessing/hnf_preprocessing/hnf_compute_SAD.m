function [rapp,var_front,var_back,mean_front,mean_back]=hnf_compute_SAD(topog,chanlocs,n)
% computeSAD() - Computes Spatial Average Difference feature
%
% Usage:
%   >> [rapp,var_front,var_back,mean_front,mean_back]=computeSAD(topog,chanlocs,n);
%
% Inputs:
%   topog      - topographies vector
%   chanlocs   - EEG.chanlocs struct
%   n          - number of ICs
%   nchannels  - number of channels
%
% Outputs:
%   rapp       - SAD values
%   var_front  - Frontal Area variance values
%   var_back   - Posterior Area variance values
%   mean_front - Frontal Area average values
%   mean_back  - Posterior Area average values
%
%

nchannels=length(chanlocs);

% Define scalp zones

% Find electrodes in Frontal Area (FA)
dimfront=0; %number of FA electrodes
index1=zeros(1,nchannels); %indexes of FA electrodes

for k=1:nchannels
    if (abs(chanlocs(1,k).theta)<60) && (chanlocs(1,k).radius>0.40) %electrodes are in FA
        dimfront=dimfront+1; %count electrodes
        index1(1,dimfront)=k;
    end
end

% Find electrodes in Posterior Area (PA)
dimback=0;
index3=zeros(1,nchannels);
for h=1:nchannels
    if (abs(chanlocs(1,h).theta)>110)
        dimback=dimback+1;
        index3(1,dimback)=h;
    end
end

if dimfront*dimback==0
    % modification Eric Rawls 9/10/2026 - return warning instead of
    % hard error
    rapp=[];
    var_front=[];
    var_back=[];
    mean_front=[];
    mean_back=[];
    warning('hnf_computeSAD:NoChannels', 'No channels included in some scalp areas. SAD computation aborted.');
    return
end

% Outputs

rapp=zeros(1,n); % SAD
mean_front=zeros(1,n); % FA electrodes mean value
mean_back=zeros(1,n); % PA electrodes mean value
var_front=zeros(1,n); % FA electrodes variance value
var_back=zeros(1,n); % PA electrodes variance value

% Output computation

for i=1:n % for each topography

    %create FA electrodes vector
    front=zeros(1,dimfront);
    for h=1:dimfront
        front(1,h)=topog(i,index1(1,h));
    end

    %create PA electrodes vector
    back=zeros(1,dimback);
    for h=1:dimback
        back(1,h)=topog(i,index3(1,h));
    end



    %compute features

    rapp(1,i)=abs(mean(front))-abs(mean(back)); % SAD
    mean_front(1,i)=mean(front);
    mean_back(1,i)=mean(back);
    % change Eric Rawls 9/10/2026 - modified sample to population variance
    % according to manuscript equation 5
    var_back(1,i)=var(back,1);
    var_front(1,i)=var(front,1);

end

