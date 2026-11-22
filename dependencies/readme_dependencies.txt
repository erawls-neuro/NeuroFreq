dependencies/ README
——————————----------------------------------------------------------------
The dependencies folder of NeuroFreq will contain external functions
that implement key operations used in NeuroFreq. These external libraries 
are generally used for one of several reasons:
1) the library has better methods for something than I can write;
2) the library is extremely well-vetted and reliable
Run nf_install from the NeuroFreq root to check these packages and the
EEGLAB plugins used by individual NF functions. It can locate an existing
installation or download a copy into dependencies/_nf_local. The package
catalog and exact function probes are in private/nf_catalog.m.

Upstream sources for the original five dependencies:

ept-tfce:      https://github.com/mensen/ept_tfce-matlab
FMINSEARCHBND: https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon
HAPPE:         https://github.com/PINE-Lab/HAPPE
MADE:          https://github.com/ChildDevLab/MADE-EEG-preprocessing-pipeline
fieldfun:      https://github.com/WD40andTape/fieldfun
