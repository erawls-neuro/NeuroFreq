function test_nf_signal_range()
% Verify that representable scientific results survive intermediate scales.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; oldRng = rng;
restore = onCleanup(@() cleanup(oldPath,oldRng)); %#ok<NASGU>
for d = {'timefrequency','spectralparameterization','phaseamplitudecoupling','setutilities'}
    addpath(genpath(fullfile(root,d{1})));
end
f = [1e-200 1 1e200]; p = [2 0 5];
assert(max(abs(hnf_expo_function(f,p)-(2-5*log10(f))))<1e-12);
p(2) = 1;
truth = [2,2-log10(2),-998];
assert(max(abs(hnf_expo_function(f,p)-truth))<1e-12);
assert(hnf_error_expo_function(p,f,truth)<1e-22);
if hnf_mex_available('hnf_knee_loss_mex')
    assert(hnf_knee_loss_mex(p,log(f),truth)<1e-22);
end
assert(abs(hnf_specparam_slope_guess([1e-200 1e200],[1e308 -1e308])/5e305-1)<1e-14);
near = [1e250,1e250+eps(1e250)];
assert(isfinite(hnf_specparam_slope_guess(near,[2 1])));
% Although the unit-height Lorentzian tail underflows, its weighted height
% is 1e-100. Recover the same SSE (and a finite gradient) with the MEX.
for kind={'gaussian','cauchy'}
    type=kind{1}; cauchy=strcmp(type,'cauchy');
    if cauchy, pts=[1e200 2e200]; expected=[1e-100 2.5e-101];
    else, pts=[40 41]; expected=exp(log(1e300)-.5*pts.^2); end
    pars=[0 1e300 1];
    actual=hnf_specparam_peak_model(pts,pars,type);
    assert(max(abs(actual./expected-1))<2e-13);
    [s,g]=hnf_error_model_constr(pars,pts,expected,type);
    assert(all(isfinite(g)));
    if hnf_mex_available('hnf_peak_objective_mex')
        [v,d]=hnf_peak_objective_mex(pars,pts,expected,[],0,double(cauchy));
        assert(all(isfinite(d)) && abs(s-v)<1e-110);
    end
    [s,g]=hnf_error_model_constr([0 1 1e-200],[0 1],[1 0],type);
    assert(isfinite(s) && all(isfinite(g)));
end

% Same slope and R-squared across ~500 orders of magnitude in log power.
f = 1:64; y = 2 - .7*log10(f);
for scale=[1e-240 1 1e240]
    fit = nf_specparam(y*scale,f,'maxPeaks',0);
    assert(abs(fit.aperiodicParms(2)/scale-.7)<2e-12 && abs(fit.RSquare-1)<1e-12);
    % Force the production fixed-spectrum MEX dispatch threshold.
    tf=struct('nsensor',1,'ntrls',1,'scale','log10','freqs',f,'power',repmat(y(:)*scale,1,1024));
    fit=nf_tfspecparam(tf,'maxPeaks',0);
    assert(max(abs(fit.specparam.ap_exponent(:)/scale-.7))<2e-12);
    assert(all(abs(fit.specparam.rsquare(:)-1)<1e-12));
end

% The output power surface is double even when coefficients are single.
Fs=128;t=(0:255)/Fs;x=cos(2*pi*8*t);
for amp=[1e-25 1e25]
    s=nf_stft(single(x*amp),Fs,1,50,1);
    assert(abs(s.power(9,1)/amp^2-1)<1e-3);
end

if hnf_mex_available('hnf_powphase_mex')
    for z={single(complex([0 3e25],[0 4e25])),complex([2e-160 -2],[3e-160 0])}
        [power,phase] = hnf_powphase_mex(z{1});
        ref = abs(double(z{1})).^2;
        assert(all(abs(power-ref)<=3*eps(ref)));
        assert(max(abs(phase-angle(z{1})))<=3*eps(max(abs(phase))));
    end
end
% Log-power 400 corresponds to an amplitude of 1e200; a direct 10^400
% intermediate used to turn this into NaN in a power-weighted estimator.
tf=struct('nsensor',1,'ntrls',1,'Fs',128,'times',t,'freqs',[4 22], ...
    'scale','log10','power',400*ones(2,numel(t)),'phase',zeros(2,numel(t)));
s=nf_pac(tf,[4 4],[22 22],'method','MVL','weighting','amplitude','nperm',0);
assert(abs(s.PAC/1e200-1)<1e-12);

% MVL permutation z scores are invariant to the units of the envelope.
rng(918); ph=rand(48,3)*2*pi-pi; weights=rand(48,3);
for cross=[false true]
    for native=[false true]
        rng(123); obs=hnf_mvl_vec(ph,weights);
        ref=hnf_pac_permutation_z(ph,weights,obs,20,'mvl',0,native,cross);
        for scale=[1e-220 1e220]
            rng(123); z=hnf_pac_permutation_z(ph,weights*scale,obs*scale,20,'mvl',0,native,cross);
            assert(max(abs(ref(:)-z(:)))<1e-12);
        end
    end
end
a = ((1:50)'*1e50); w=mod((1:50)',7)+1;
direction=atan2(sin(a),cos(a));
assert(abs(hnf_mi_vec(a,w,4)-hnf_mi_vec(direction,w,4))<1e-12);
if hnf_mex_available('hnf_mi_mex')
    assert(abs(hnf_mi_mex(a,w,4)-hnf_mi_vec(direction,w,4))<1e-12);
end
fprintf('Spectral dynamic range, peak objectives, PAC scaling and TF power passed.\n');
end

function cleanup(p,r)
path(p);rng(r);
end
