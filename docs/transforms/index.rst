
transforms
==========

The primary functionality of NeuroFreq is to perform efficient, high-quality TF decompositions with a wide variety of algorithms. Following optional data preparation, M/EEG time series data can be transformed to a TF representation using any of the included algorithms, including several implementations of linear and quadratic decomposition methods.

Different TF transformation algorithms are characterized by their time-frequency tradeoff. STFT has an inflexible TF tradeoff that can favor either time or frequency resolution but not both. With CWT, DCWT, filter-Hilbert, and S-transform, there is good temporal resolution at high frequencies and good frequency resolution at lower frequencies. The Wigner-Ville distribution, a prototypical TF distribution, has excellent joint time-frequency resolution, but suffers from cross-terms that render it unusable for most practical applications. As such it is not included in NeuroFreq. The Reduced Interference Distributions solve this problem by applying kernel functions to the Wigner-Ville distribution, significantly masking cross-terms at the expense of slightly lower TF resolution.

.. image:: Fig_TFresolution.png
  :width: 600

Each transform in NeuroFreq is demonstrated on a set of synthetic data made publicly available by Arts & van den Broke (2022). We appreciate the authors making these synthetic data available to the community.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   nf_stft/index
   nf_filterhilbert/index
   nf_demodulation/index
   nf_dcwt/index
   nf_cwt/index
   nf_stransform/index
   nf_ridbinomial/index
   nf_ridbornjordan/index
   nf_ridrihaczek/index


References
^^^^^^^^^^
Arts, L.P. and van den Broek, E.L., 2022. The fast continuous wavelet transformation (fCWT) for real-time, high-quality, noise-resistant time–frequency analysis. Nature Computational Science, 2(1), pp.47-58.