NeuroFreq/ README
——————————----------------------------------------------------------------
The preproc folder of NeuroFreq contains functions to implement the most 
common preprocessing methods in EEG/neurophysiology research. These 
functions include several preset pipelines (PREP, FASTER, cleanrawdata, 
MADE, & our in-house BDC), as well as fully customizable filtering, bad 
channel detection, ICA cleaning routines, pre cleaning via GEDAI or ASR,
and final thresholding and optional epoch-level interpolation. 

External pipeline presets require installation of the pipeline requested 
as we do not re-implement existing pipeline code. You can find the
required pipeline scripts at: 

PREP:          https://github.com/VisLab/EEG-Clean-Tools
FASTER:        https://sourceforge.net/projects/faster/files/
clean_rawdata: https://github.com/sccn/clean_rawdata
HAPPE:         https://github.com/PINE-Lab/HAPPE

Functions:

1) nf_badchans:
Configurable bad channel detection and deletion. 

2) nf_cleanic: 
Configurable ICA-based artifact removal and cleaning.

3) nf_correctlatencies:
Correct event onset latencies using another set of events. Intended for
Event correction paired with a device like Cedrus StimTracker, Cortech 
Triggy, Black Box Toolkit, or similar.

4) nf_eegquality:
Computes a pre-post artifacting QC dashboard for an EEG dataset upon 
request during preprocessing (not to be called independently).

5) nf_epochbynumber:
A simplified method of using EEGLAB's epoch function to select only 
certain epochs, originally written by Adrian Fischer.

6) nf_filter:
Configurable filtering and downsampling of raw/continuous EEG data.

7) nf_preprocess:
A primary wrapper function that completely preprocesses an input EEG
dataset using various presets and configurations. Calls filter, 
badchans, cleanic, and thresh in order. 

8) nf_readpsychopy:
Reads a psychopy behavioral file using some heuristics to determine which
rows contain real data.

9) nf_thresh:
Configurable thresholding and deletion/interpolation to complete 
preprocessing. 