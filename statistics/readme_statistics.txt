statistics/ README
——————————----------------------------------------------------------------
The stats folder of NeuroFreq contains functions for running statistics
with time-frequency data structures.

Functions:

1) nf_corr_TFCE
Uses TFCE to correlate EEG with an external variable. Can be single trials
or averaged subject level data and variables, as long as dimensions are
correct. For multivariate correlations (e.g. corr(subjectXchannelXtime,
subjectXchannelXtime), uses very efficient vectorization techniques.
returns a standard 'stats' structure that can be plotted.

2) nf_glm_TFCE
Uses TFCE to regress EEG onto an set of external variables. Can be single 
trials or averaged subject level data and variables, as long as dimensions 
are correct. returns a standard 'stats' structure that can be plotted.

3) nf_indsamples_TFCE
Uses TFCE to t-test EEG against another EEG set. Can be single trials
or averaged subject level data, as long as dimensions are
correct. returns a standard 'stats' structure that can be plotted.

4) nf_lme_TFCE
Uses TFCE to correct results of a mixed-effects regression analysis. Must
be aggregated single-trial data across multiple subjects, and the Subject
ID variable must be identified. A random intercept per subject is enforced
and random slopes can be added per regressor. Returns a standard 'stats'
structure that can be plotted.

5) nf_makefactorialdesign
Makes a factorial design matrix from TF.behavior. Includes interaction 
terms as specified and optional subject ID for LME.

6) nf_onesample_TFCE
Uses TFCE to t-test EEG against h0. Can be single trials
or averaged subject level data, as long as dimensions are
correct. returns a standard 'stats' structure that can be plotted.

