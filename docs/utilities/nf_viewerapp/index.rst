
nf_viewerapp
============

NeuroFreq includes a viewer app that allows interacting with TF structures for M/EEG including topographic plots, surface plots, and scrolling through single trials. Since the TF structure requires channel locations for plotting, this function only applies to TF computed using the nf_tftransform function. The function builds on the excellent tfViewerX.m function from Mike X. Cohen's 2014 book Analyzing Neural Time Series. Additions to the function include scrolling through multiple trials/conditions, and line plots of time/frequency average activities. The viewing app is called as follows:

.. code-block:: matlab
   
  nf_viewerapp(TF);

Where 1) TF is the output of nf_tftransform (single subject) or of nf_aggregate (multi-subject).

References
^^^^^^^^^^
Cohen, M. X. (2014). Analyzing neural time series data: theory and practice. MIT press.