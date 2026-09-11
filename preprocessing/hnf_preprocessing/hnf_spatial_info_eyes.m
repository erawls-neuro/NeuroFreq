function ICs_to_keep = hnf_spatial_info_eyes(blinks,eye_electrodes,topog,chanlocs,n)
% HNF_SPATIAL_INFO_EYES   Calculate spatial info on eye movements for developmental data
%
% Note: this function was rewritten in large part by Eric Rawls 9/2026.
% Check intended use case before applying.
%

%define variables
nchannels=length(chanlocs);
ICs_to_keep = [];
arb_dist = 5.5;
arb_zval = 2;

% Eric Rawls modified 9/10/2026.
% Standardize each component across all channels; allow non-square ICA.
zmatrix = zeros(n,nchannels);
for i = 1:n
    curr_vect = topog(i,1:nchannels).';
    z_vect = zscore(curr_vect);
    zmatrix(i,:) = z_vect.';
end

%list of channels on the outside of the head
% Find electrodes in outer ring
dimouterring=0; %number of electrodes
indexor=zeros(1,nchannels); %indices of electrodes
for k=1:nchannels
    if (chanlocs(1,k).radius>0.51)
        dimouterring=dimouterring+1; %count electrodes
        indexor(1,dimouterring)=k;
    end
end
outer_ring = nonzeros(indexor);
outer_ring = unique([outer_ring(:); eye_electrodes(:)]).';
dimhigherthresh=0;
indexht=zeros(1,nchannels);
for k=1:nchannels
    if (chanlocs(1,k).radius>0.35 && chanlocs(1,k).radius<0.5 && abs(chanlocs(1,k).theta)<20)
        dimhigherthresh=dimhigherthresh+1;
        indexht(1,dimhigherthresh)=k;
    end
end
higherthresh = nonzeros(indexht)';

keep_IC = []; % will be the list of ICs to keep
for ic=1:n
    if sum(ic==blinks)>0 % if this IC matches any numbers in the eye artifact list
        % find any channels/electrodes that have high activity (z>2)
        chans=find(abs(zmatrix(ic,:))>arb_zval); %grabs any chans/electrodes that are above thresh
        % delete the outer ring(s) of electrodes
        % We only care about the spread over more central sites for eye stuff
        chans_to_delete=[];
        for e=1:length(chans)
            if find(outer_ring == chans(e)) > 0
                % delete this chan/electrode from the list of chans (we don't care about it)
                chans_to_delete=[chans_to_delete e]; %#ok <AGROW>
            end
            if find(higherthresh == chans(e)) > 0
                % Eric Rawls modified 9/10/2026: restore ICA sign invariance.
                % The original signed comparison discarded negative channels even
                % when their magnitude exceeded 2.5 (e.g., -3 versus +3).
                if abs(zmatrix(ic,chans(e))) < 2.5
                    % this channel isn't above the higher threshold set for this ring on the net
                    chans_to_delete=[chans_to_delete e]; %#ok <AGROW>
                end
            end
        end
        chans(unique(chans_to_delete))=[];

        % calculate the distance between the points
        pts_matrix_3D=[];
        for c=1:length(chans) % loop through list of chans above the threshold value and get a list of points
            pts_matrix_3D=[pts_matrix_3D; chanlocs(1,chans(c)).X chanlocs(1,chans(c)).Y chanlocs(1,chans(c)).Z]; %#ok <AGROW>
        end
        dist_matrix = pdist(pts_matrix_3D,'euclidean'); % calculate distances btw all pts. above the threshold

        % figure out if component should be kept
        corresponding_pts = []; % gives us the indices of the points
        for jj=1:length(chans)-1
            corresponding_pts = [corresponding_pts; (jj*ones(1,length(chans)-jj))' ((jj+1):length(chans))']; %#ok <AGROW>
        end
        dists_below_thresh = dist_matrix<arb_dist; % tells us distances below the threshold
        pts_within_dist = corresponding_pts(dists_below_thresh,:); % tells us which point pairs are within the distnace threshold
        if ~isempty(pts_within_dist) % if any pts are within the arbitrary distance
            most_com = mode(pts_within_dist(:)); % lets us know if eye stuff spreads too far into the head
            if length(find(pts_within_dist(:) == most_com))>1
                keep_IC = [keep_IC ic];  %#ok <AGROW>
            end
        end
        ICs_to_keep = unique(keep_IC);
    end
end

end