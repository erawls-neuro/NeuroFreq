%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HNF_EM - ADJUST EM with minimal NeuroFreq corrections
% 
% Performs automatic threshold on the digital numbers 
% of the input vector 'vec'; based on Expectation - Maximization algorithm
%
% Eric Rawls modified 9/11/2026 for NeuroFreq internal testing.
% This file starts from the original ADJUST EM function. The original
% initialization, class updates, costs, and section layout are retained.
% Changes below address numerical stability, convergence, and the agreed
% warning/continuation behavior. The variance floor and root-selection
% correction can change classifications; they are described where applied.

% Reference paper:
% Bruzzone, L., Prieto, D.F., 2000. Automatic analysis of the difference image 
% for unsupervised change detection. 
% IEEE Trans. Geosci. Remote Sensing 38, 1171:1182

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage:
%   >> [last,med1,med2,var1,var2,prior1,prior2]=hnf_EM(vec);
%
% Input: vec (finite real row or column vector, to be thresholded)
%
% Outputs: last (threshold value)
%          med1,med2 (mean values of the Gaussian-distributed classes 1,2)
%          var1,var2 (variance of the Gaussian-distributed classes 1,2)
%          prior1,prior2 (prior probabilities of the Gaussian-distributed classes 1,2)
% Class 1 is returned with the higher mean. If no usable fit/threshold
% remains, warn and return NaN for all outputs. ADJUST's comparisons then
% skip only the decisions requiring this threshold. No caller changes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2009-2014 Andrea Mognon (1) and Marco Buiatti (2), 
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


function [last,med1,med2,var1,var2,prior1,prior2]=hnf_EM(vec)

% Eric Rawls modified 9/11/2026: reject malformed inputs, but handle expected
% fitting failures locally below. Do not silently drop IC feature values.
if ~isnumeric(vec) || ~isreal(vec) || ~isvector(vec) || numel(vec)<2 ...
        || any(~isfinite(vec(:)))
    error('EM:InvalidInput','Input must be a finite real vector with at least two values.');
end
vec=double(vec);
featureName=inputname(1);
if isempty(featureName), featureName='input feature'; end
try

if size(vec,2)>1
	len=size(vec,2); %number of elements
else
	vec=vec';
	len=size(vec,2); 
end

c_FA=1; % False Alarm cost
c_MA=1; % Missed Alarm cost

med=mean(vec);
standard=std(vec);
% Same midrange as ADJUST; divide first to avoid overflow in the sum.
mediana=max(vec)/2+min(vec)/2;

alpha1=0.01*(max(vec)-mediana); % initialization parameter/ righthand side
alpha2=0.01*(mediana-min(vec)); % initialization parameter/ lefthand side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPECTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

train1=[]; % Expectation of class 1
train2=[];
train=[]; % Expectation of 'unlabeled' samples

for i=1:(len)
    if (vec(i)<(mediana-alpha2)) 
        train2=[train2 vec(i)];
    elseif (vec(i)>(mediana+alpha1))
        train1=[train1 vec(i)];
    else
	  train=[train vec(i)];
    end
end

n1=length(train1);
n2=length(train2);
if n1==0 || n2==0
    error('EM:DegenerateFit','The initialization does not define two nonempty classes.');
end

med1=mean(train1);
med2=mean(train2);
prior1=n1/(n1+n2);
prior2=n2/(n1+n2);
var1=var(train1);
var2=var(train2);

% Eric Rawls modified 9/11/2026: the original Bayes helper returned density
% 1 when variance was zero. Use a small positive variance instead, scaled
% to this feature's variance. Apply the same floor after every update.
% 1e-6 is our explicit regularization choice, not a published ADJUST value.
varianceScale=var(vec,1);
if ~isfinite(varianceScale) || varianceScale<=0
    error('EM:NumericalFailure','Feature variance is not finite and positive.');
end
varianceFloor=1e-6*varianceScale;
regularized=any([var1 var2]<varianceFloor);
if var1<varianceFloor, var1=varianceFloor; end
if var2<varianceFloor, var2=varianceFloor; end
initial=[med1 med2 var1 var2 prior1 prior2];
if any(~isfinite(initial)) || min([var1 var2])<=0
    error('EM:DegenerateFit','Both seed classes must have finite, positive variance.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAXIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=0;
dif_med_1=Inf; % difference between current and previous mean
dif_med_2=Inf;
dif_var_1=Inf; % difference between current and previous variance
dif_var_2=Inf;
dif_prior_1=Inf; % difference between current and previous prior
dif_prior_2=Inf;
% Eric Rawls modified 9/11/2026: scale tolerances by the initial statistics,
% as in the reviewed formulation. The roundoff floor handles zero means.
% Inf above ensures that the first iteration runs even for large tolerances.
stop=max(0.0001*abs(initial), ...
    32*eps*[max(abs(vec)) max(abs(vec)) varianceScale varianceScale 1 1]);
max_iter=1000;

% Original && stopped when ANY parameter stabilized. Continue while ANY
% change is too large, so all six must stabilize; cap the iteration count.
while any([dif_med_1 dif_med_2 dif_var_1 dif_var_2 ...
        dif_prior_1 dif_prior_2]>stop) && count<max_iter

    count=count+1;

    med1_old=med1;
    med2_old=med2;
    var1_old=var1;
    var2_old=var2;
    prior1_old=prior1;
    prior2_old=prior2;
	prior1_i=zeros(1,len);
	prior2_i=zeros(1,len);

    % FOLLOWING FORMULATION IS ACCORDING TO REFERENCE PAPER:
    
    for i=1:len
        % Eric Rawls modified 9/11/2026: same responsibility calculation in
        % log space. Subtract the larger log density before exponentiating
        % so ordinary density underflow cannot produce 0/0.
        log1=log(prior1_old)+logBayes(med1_old,var1_old,vec(i));
        log2=log(prior2_old)+logBayes(med2_old,var2_old,vec(i));
        logmax=max(log1,log2);
        if ~isfinite(logmax)
            error('EM:NumericalFailure','Gaussian log densities exceed floating-point range.');
        end
        p1=exp(log1-logmax);
        p2=exp(log2-logmax);
        prior1_i(i)=p1/(p1+p2);
        prior2_i(i)=p2/(p1+p2);
    end
	
	
	prior1=sum(prior1_i)/len;
	prior2=sum(prior2_i)/len;
    if any(~isfinite([prior1 prior2])) || min([prior1 prior2])<=0
        error('EM:DegenerateFit','A fitted class has zero or nonfinite responsibility mass.');
    end
	med1=sum(prior1_i.*vec)/(prior1*len);
	med2=sum(prior2_i.*vec)/(prior2*len);
    % Retain the original/published previous-mean variance update.
	var1=sum(prior1_i.*((vec-med1_old).^2))/(prior1*len);
	var2=sum(prior2_i.*((vec-med2_old).^2))/(prior2*len);
    regularized=regularized || any([var1 var2]<varianceFloor);
    if var1<varianceFloor, var1=varianceFloor; end
    if var2<varianceFloor, var2=varianceFloor; end
    if any(~isfinite([med1 med2 var1 var2])) || min([var1 var2])<=0
        error('EM:DegenerateFit','A fitted class has nonfinite parameters or nonpositive variance.');
    end

    dif_med_1=abs(med1-med1_old);
    dif_med_2=abs(med2-med2_old);
    dif_var_1=abs(var1-var1_old);
    dif_var_2=abs(var2-var2_old);
    dif_prior_1=abs(prior1-prior1_old);
    dif_prior_2=abs(prior2-prior2_old);

end

% Eric Rawls modified 9/11/2026: report regularization once per fit. At the
% iteration limit, retain valid final estimates and warn instead of stopping
% the dataset. This continues processing without claiming convergence.
if regularized
    warning('EM:VarianceRegularized', ...
        '%s: applied variance floor 1e-6 times the whole-feature variance.',featureName);
end
if any([dif_med_1 dif_med_2 dif_var_1 dif_var_2 dif_prior_1 dif_prior_2]>stop)
    warning('EM:NoConvergence', ...
        ['%s: not all EM parameters stabilized within %d iterations. ' ...
         'Following legacy ADJUST''s permissive behavior, continuing with ' ...
         'the final parameter estimates.'],featureName,max_iter);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THRESHOLDING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eric Rawls modified 9/11/2026: keep the means, variances, and probabilities
% paired if the class order reverses. Class 1 must be the higher-mean class.
if med1<med2
    [med1,med2]=deal(med2,med1);
    [var1,var2]=deal(var2,var1);
    [prior1,prior2]=deal(prior2,prior1);
end
gap=med1-med2;
if ~isfinite(gap) || gap<=0
    error('EM:NoSeparatingThreshold','The fitted means do not define distinct ordered classes.');
end

% Same ADJUST density-intersection equation, centered/scaled only during
% root calculation to avoid subtracting large squared means. The fitted
% parameters and returned threshold retain their original units.
center=med2+gap/2;
scale=max([sqrt(var1) sqrt(var2) gap]);
u1=(med1-center)/scale;
u2=(med2-center)/scale;
v1=(sqrt(var1)/scale)^2;
v2=(sqrt(var2)/scale)^2;
if min([v1 v2])<=0
    error('EM:NumericalFailure','Variance ratios exceed floating-point range.');
end
k=c_MA/c_FA;
logratio=log(k)+log(prior1)-log(prior2)+0.5*(log(v2)-log(v1));
a=(v1-v2)/2;
b=v2*u1-v1*u2;
c=v1*v2*logratio+(v1*u2^2-v2*u1^2)/2;
coeffScale=max(abs([a b c]));
if ~isfinite(coeffScale) || coeffScale==0
    error('EM:NoSeparatingThreshold','No numerically identifiable density intersection.');
end
a=a/coeffScale;
b=b/coeffScale;
c=c/coeffScale;

% Eric Rawls modified 9/11/2026: equal variances give a linear equation.
% For unequal variances, calculate the original plus-sign root without
% cancellation. The factored discriminant also avoids loss of precision.
if a==0
    if b==0
        error('EM:NoSeparatingThreshold','The density equation has no unique crossing.');
    end
    soglia1=-c/b;
else
    rad=(u1-u2)^2-2*(v1-v2)*logratio;
    if rad<=0
        error('EM:NoSeparatingThreshold','The weighted densities have no distinct real crossing.');
    end
    sqrtRad=(sqrt(v1)*sqrt(v2)/coeffScale)*sqrt(rad);
    s=1;
    if b<0, s=-1; end
    q=-0.5*(b+s*sqrtRad);
    if b<0
        soglia1=q/a;
    else
        soglia1=c/q;
    end
end

% The original selected the other root when soglia1 was outside the means.
% That other root crosses toward the lower-mean class. Keep the rising root
% for ADJUST's feature>=threshold rule, even outside the means, and warn.
last=center+scale*soglia1;
if ~isfinite(last)
    error('EM:NumericalFailure','The density crossing is not finite.');
end
if soglia1<u2 || soglia1>u1
    warning('EM:OutsideMeans', ...
        ['%s: no crossing lies between the fitted class means. ' ...
         'Continuing with the crossing toward the higher-mean class.'],featureName);
end

catch failure
    % Eric Rawls modified 9/11/2026: the old NaN fallback used the feature
    % midrange, and a negative discriminant could leave last unassigned.
    % If a usable threshold still cannot be calculated, warn and return NaN
    % so independent ADJUST detectors can continue. Rethrow actual code errors.
    if ~any(strcmp(failure.identifier, ...
            {'EM:DegenerateFit','EM:NumericalFailure','EM:NoSeparatingThreshold'}))
        rethrow(failure);
    end
    [last,med1,med2,var1,var2,prior1,prior2]=deal(NaN);
    warning(failure.identifier, ...
        '%s: %s Returning NaN; dependent ADJUST decisions are unavailable.', ...
        featureName,failure.message);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prob=logBayes(med,var,point)
% Eric Rawls modified 9/11/2026: return log density for stable responsibilities.
% The common -log(2*pi)/2 cancels. The caller now floors zero/tiny variances.
prob=-0.5*(log(var)+((point-med)/sqrt(var))^2);
end
