function beta = hnf_lsreg_pinv(D,regs)
%
% GENERAL
% -------
% calculate least-squares regression using Moore-Penrose pseudoinverse.
% dim 1 of data/regs should be subjects/observations, dim2 should be 
% variables 
%
% INPUTS
% ------
% D - matrix of data (DVs). n(subjects) X n(variables).
% regs - matrix of regressors (IVs). n(subjects) X n(variables)
%
% OUTPUT
% ------
% beta - matrix of beta weights resulting from least-squares fit
%
% -----
% E. Rawls, erawls89@gmail.com, rawls017@umn.edu. 
% Nov. 2023
% Copyright (c) 2023 by E. Rawls.
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of version 3 of the GNU General Public License as 
% published by the Free Software Foundation.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

%standardize
D = (D - mean(D)) ./ std(D);
regz = (regs - mean(regs)) ./ std(regs);
%make column of 1s for regression
regz = [ones(size(regz,1),1), regz];
%fit glm using pseudoinverse
beta = pinv(regz) * D; %moore-penrose pseudo-inverse
%drop intercept
beta = beta(2:end,:);
end