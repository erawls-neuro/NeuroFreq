function test_nf_octave_tfce()
% Known midpoint integrals for connected/disconnected, mixed-sign maps.
root = fileparts(fileparts(mfilename('fullpath')));
p = path;
addpath(fullfile(root,'statistics','hnf_statistics'), ...
    fullfile(root,'phaseamplitudecoupling','hnf_phaseamplitudecoupling'));
state = rng;
policy = hnf_mex_available('scope',struct());
restore = onCleanup(@() restore_state(p,state,policy)); %#ok<NASGU>
% EPT's positive traversal used to overwrite the next negative seed.
signed = hnf_tfce_transformation([1 -1],1,[1 0],false);
assert(max(abs(signed-[1 -1]))<1e-14);
maps = {[2 0 0;0 2 -2], cat(3,[2 0;0 0],[0 -2;0 2])};
for i=1:numel(maps)
    x = maps{i};
    integral = sum((((0:49)+.5)*2/50).^2)*2/50;
    expected = (double(x>0)*2^.66-double(x<0))*integral;
    enhanced = hnf_tfce_transformation(x,[1 2;2 1],[.66 2],false);
    assert(max(abs(enhanced(:)-expected(:)))<1e-12);
    reversed = hnf_tfce_transformation(-x,[1 2;2 1],[.66 2],false);
    assert(max(abs(reversed(:)+enhanced(:)))<1e-12);
end
% Adjacent points split at h=1. Isolated rows stay separate. A zero
% terminates the row even when a nonzero follows it.
x = [2 1;0 -2];
y = hnf_tfce_integrate(x,[1 0 2;2 0 1],[1 0]);
assert(max(abs(y(:)-[3;0;2;-2]))<1e-12);
% Directed rows use the explicit traversal, following column-major seeds.
y = hnf_tfce_integrate([1;1],[1 2;2 0],[1 0]);
assert(max(abs(y(:)-2))<1e-12);
y = hnf_tfce_integrate([1;1],[1 0;2 1],[1 0]);
assert(max(abs(y(:)-1))<1e-12);
% Exact midpoint extent 8 times the squared-height integral gives
% 2.6664e-321, or 540 units of the least double subnormal. A large
% extent power may itself overflow even when the final integral fits.
check_extremes(@hnf_tfce_integrate);
% Tiny H=0 statistics still have representable signed integrals.
tiny = realmin * eps;
assert(isequal(hnf_tfce_integrate([tiny -tiny],1,[0 0]), [tiny -tiny]));
% Dense symmetric, directed, repeated, zero-terminated, and empty rows.
if exist('hnf_tfce_clusters_mex','file') == 3
    check_extremes(@hnf_tfce_clusters_mex);
    rng(114,'twister');
    for rows = {[1 2 2;2 1 0;3 0 3],[0 2 1;1 2 0;2 3 0]}
        for sz = {[3 8],[3 4 7],[3 1 4]}
            x = round(randn(sz{1})*8)/8;
            for e = {[.66 2],[1 0],[0 0]}
                oracle = hnf_tfce_integrate(x,rows{1},e{1});
                hnf_mex_available('scope',[]);
                value = hnf_tfce_integrate(x,rows{1},e{1});
                hnf_mex_available('scope',struct());
                assert(isequal(size(oracle),size(value)));
                assert(all(abs(oracle(:)-value(:)) <= 2e-13*(1+abs(oracle(:)))));
            end
        end
    end
    % Guard the neighbouring and dimension-index ABIs and seed sign.
    x = [realmax/100 -realmax/100 0];
    assert(max(abs(hnf_tfce_clusters_mex(x,1,[0 0])-x)) <= max(abs(x))*3e-15);
    assert(isequal(hnf_tfce_clusters_mex([tiny -tiny],1,[0 0]), [tiny -tiny]));
    assert(isequal(hnf_tfce_clusters_mex(zeros(0,2),zeros(0,1),[0 0]),zeros(0,2)));
    x = zeros(1,40000); x(end-1:end) = [1 -1];
    assert(max(abs(hnf_tfce_clusters_mex(x,1,[0 0])-x))<1e-14);
    bad = {@() hnf_tfce_clusters_mex([1 2],2,[1 0]), ...
           @() hnf_tfce_clusters_mex(Inf,1,[0 0]), ...
           @() hnf_tfce_clusters_mex([1 2],1,[0 -1])};
    for k=1:numel(bad)
        caught=false;
        try, unused=bad{k}(); catch err, caught=strncmp(err.identifier,'NeuroFreq:',10); end %#ok<NASGU>
        assert(caught);
    end
end
fprintf('TFCE signed 2D/3D integrals and adjacency checks passed.\n');
end

function restore_state(p,state,policy)
hnf_mex_available('scope',policy); rng(state); path(p);
end

function check_extremes(fn)
v = fn(ones(1,8)*1e-107,1,[1 2]);
assert(all(v/(realmin*eps) == 540));
v = fn(ones(1,2)*2^-1000,1,[2000 0]);
assert(all(abs(v/2^1000-1)<2e-12));
end
