utilities/ README
——————————----------------------------------------------------------------
The utilities folder of NeuroFreq contains high-level functions for working 
directly with EEG and time-frequency data structures. These functions 
include utilities to prepare data for TF analysis, run TF analysis, and 
manipulate TF surfaces including resampling, averaging, and baseline 
correction.

Functions:

1) nf_aggregate: 
Opens a directory window to select TF sets. Aggregates multiple data files 
into a single TF structure for group analysis.

2) nf_avebase:
Averages and optionally baseline corrects single trial TF structures. 
Trials are averaged together per sensor, and an optional dB/z-score/percent 
baseline correction is applied. Phase data are averaged as the inter-trial
phase coherence (ITPC) measuring phase consistency over trials.

3) nf_cpm:
Implements the concurrent phaser method for separating induced and evoked
Time frequency power. Averages trials as part of the CPM definition. Trials
are averaged into phased-locked and non-phase-locked activities.

4) nf_makebehavior
Creates a standard NeuroFreq behavior structure from a csv file of a task
(for example, stimuli, RTs, accuracy, etc.). Must have column headers. Adds
.behavior field to TF structure which is automatically carried forward
through further processing.

5) nf_prepdata:
Accepts either an input EEGLAB .set in memory or a 1/2/3D data tensor. In 
either case, the function will 1) remove quadratic trends 
from single-trial data, 2) center the data, 3) cosine-square taper 
single-trial segments, and 4) convert signal to analytic using the Hilbert 
transform. If an EEGLAB .set is input, then a prepared EEGLAB 
.set is returned; if a data tensor is input, a prepared data tensor is 
returned.

6) nf_resample:
Resample time and frequency axes for computed TF structures. Useful for 
downsampling full-resolution surfaces for computational efficiency.

7) nf_rmerp
Removes phase-locked or evoked activity from each frequency. This
removal gets rid of any phase-locking in the data and removes contribution
of the ERP or evoked potential to TF results.

8) nf_specparam
Parameterizes power spectra into periodic and aperiodic components. 

9) nf_stregress
Computes single-trial correlation/regression on input EEG/TF data.

10) nf_tfspecparam
Uses the SPRiNT method (Wilson et al., 2022, eLife) to parameterize TF 
surfaces computed by NeuroFreq.

11) nf_tftransform:
Implements every TF transform from the tf_fun folder in one 
wrapper function. Accepts only EEGLAB formatted .set files as input 
- if you want to analyze data tensors directly, use the transform functions 
directly. Read the help for full details on available keyword-argument
pairs.

12) nf_unpack2long
Unpacks a TF set to long format. Columns include 'channel', 'frequency',
'times', etc so it is highly redundant. Still the best way to extract data
for further decimation/processing/averaging in R.



