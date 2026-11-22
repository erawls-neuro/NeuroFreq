Demo: EEGLAB sample dataset
===========================

The repository includes EEGLAB's visual experiment sample and a synthetic
signal with clean and noise-corrupted versions. After ``nf_install``, run
either demo from MATLAB:

.. code-block:: matlab

  root = fileparts(which('nf_install'));
  addpath(fullfile(root,'demo','timefrequency'));
  demo_EEGLAB
  demo_synthetic

The EEGLAB demo removes the EOG channels, prepares the trials and compares
nine transforms. The scalp maps show mean 8–12 Hz power, while the adjacent
panels show the surface at the channel labeled Pz. For transforms with
linear oscillator phase, the second figure shows inter-trial phase
coherence (ITPC). The Rihaczek distribution has a different phase
interpretation. The synthetic demo compares both signals at the full time
scale and in a low-frequency detail view. Its quadratic transforms use the
same five-second excerpt to limit memory use. CWT requires Wavelet Toolbox.

For a single linear transform, this function body runs within
``nf_run(@my_analysis)``:

.. code-block:: matlab

  root = fileparts(which('nf_install'));
  EEG = pop_loadset('filename',fullfile(root,'demo','timefrequency','eeglab_sample.set'));
  EEG = pop_select(EEG,'nochannel',{'EOG1','EOG2'});
  EEG = nf_prepdata(EEG);
  TF = nf_tftransform(EEG,'method','dcwt','freqs',1:30,'times',[-0.5 1]);
  avg = nf_normbase(TF,'BaselineMethod','none');
  nf_viewerapp(avg);

.. image:: nf_gui_screenshot.png
  :width: 800
