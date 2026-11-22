function [rejected, classification] = nfADJUST(EEG, varargin)
% NFADJUST  Montage-adaptive NeuroFreq implementation of ADJUST.
%
% [REJECTED, CLASSIFICATION] = NFADJUST(EEG, ...)
%
% EEG must contain an EEGLAB ICA decomposition. Unlike the released
% ADJUST implementation, the decomposition may be rectangular and may use
% a subset of EEG channels through EEG.icachansind. Spatial roles are
% transferred to the actual ICA-map montage from a normalized anatomical
% atlas; no channel labels or montage-specific indices are hard coded.
% Channel coordinates must use EEGLAB's canonical X-nose/Y-left/Z frame or
% EEGLAB theta/radius convention; nf_cleanic normalizes them automatically.
%
% Name/value inputs:
%   epochLength          continuous-data temporal-feature chunk length,
%                        seconds (default 5); existing 3-D trials remain
%                        the ADJUST temporal epochs
%   periodicRescue       rescue components with 5-15-Hz periodic structure
%                        (default false for nfADJUST)
%   spatialRescue        apply the adjusted-ADJUST spread rescue (false)
%   rescueBand           periodic rescue band in Hz (default [5 15])
%   spectrumRange        specparam fit range in Hz (default [1 40])
%   narrowPeakWidth      narrow-pass bandwidth limits (default [0.5 12])
%   plateauPeakWidth     broad-pass bandwidth limits (default [8 24])
%   emitWarnings         emit partial-evaluation warnings (default true)
%
% The returned structure separates raw artifact candidates, rescued
% components, and the final rejection set. The function has no plotting,
% report-file, current-directory, path, or random-number side effects.
%
% This implementation preserves the feature conjunctions described by
% Mognon et al. (2011) while repairing dimensional, geometry, thresholding,
% trimming, and focality calculations in the released source.
%
% SPDX-License-Identifier: GPL-2.0-or-later

[rejected, classification] = nf_adjust_core(EEG, 'adult', varargin{:});
end
