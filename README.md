# NeuroFreq
NeuroFreq is an open-source MATLAB toolbox for flexible time-frequency analysis of human brain electrophysiological fields recorded using EEG, MEG, and ECoG.

<img src="./nf_logo.png" width="200" height="200">

Time-frequency (TF) analysis is a vital tool for interpreting brain electrophysiological data including electroencephalography (EEG) and magnetoencephalography (MEG), but insights offered by TF analysis are tempered by practical challenges and inconsistencies in the application of these methods. For one example, there are many different algorithms for time-frequency analysis, and each of them has distinct strengths and weaknesses. Existing software packages largely remove the freedom of choice to go between these different transforms because they only include one or at most a few options for TF transformation. NeuroFreq includes a greater variety of time-frequency transformation algorithms than existing toolboxes, and is flexible with the capability to run both linear and quadratic time-frequency transformations.

<img src="./nf_table.png" width="500" height="200">

The NeuroFreq Toolbox is described in the accompanying preprint: https://www.biorxiv.org/content/10.1101/2023.11.01.565154v2. Read the NeuroFreq documentation at: https://neurofreq-public.readthedocs.io/en/latest/index.html. The documentation includes tutorials for common analyses as well as more information on the core functions of NeuroFreq.

## Install NeuroFreq in MATLAB

Open the NeuroFreq folder in the MATLAB desktop and run:

```matlab
report = nf_install;
```

The installer compares the files in this folder with the latest commit on NeuroFreq's GitHub default branch. If you approve an update, it verifies the download first and moves the entire old folder to a sibling named `NeuroFreq_deprecated` (numbered if necessary). Your edits, data, downloaded dependencies, and any Git history stay in that backup. The replacement is a source distribution downloaded from GitHub. If you maintain a development checkout, choose **Keep this copy** to continue working in it.

The Command Window shows which features use each package, which copy was found, and where any download went. You can point to an existing installation, keep a dependency copy with NeuroFreq, or skip a feature you do not use. For an older EEGLAB release, you can keep it, launch its updater, or install an NF-only release. Packages managed by the installer live in `dependencies/_nf_local/`. A network failure or skipped package appears in the report so you can rerun setup later.

Only the top-level NeuroFreq folder needs to stay on the MATLAB path. To run one function with the dependencies you selected, or to run an existing analysis function that contains ordinary `nf_*` calls, use:

```matlab
EEG = nf_run('nf_filter', EEG);
result = nf_run(@my_analysis, EEG);
```

`nf_run` restores the previous MATLAB path when the call returns or throws an error. To add another package to the installer, edit the catalog in [private/nf_catalog.m](private/nf_catalog.m).

# How to cite NeuroFreq
NeuroFreq is currently a preprint, and any and all input is welcome. Until NeuroFreq is published as a peer-reviewed item, please cite the preprint if you have used our software for your data analysis.