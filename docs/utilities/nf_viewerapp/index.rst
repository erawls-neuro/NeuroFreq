
nf_viewerapp
============

NeuroFreq includes a viewer app that allows interacting with TF structures for M/EEG including topographic plots, surface plots, and scrolling through single trials. Since the TF structure requires channel locations for plotting, this function only applies to TF computed using the nf_tftransform function. The function builds on the excellent tfViewerX.m function from Mike X. Cohen's 2014 book Analyzing Neural Time Series. Additions to the function include scrolling through multiple trials/conditions, and line plots of time/frequency average activities. The viewing app is called as follows:

.. code-block:: matlab
   
  viewer = nf_run('nf_viewerapp',TF);

Where 1) TF is the output of nf_tftransform (single subject) or of nf_aggregate (multi-subject).

In Octave this call opens a regular figure and returns its handle. Supply a
measure such as ``'power'``, ``'phase'`` or ``'stat'`` as the next argument.
Click the TF map to select a bin or a plotted sensor to select its map. The
arrow keys move in time and frequency. A/D move between sensors, Z/X between
observations, and M/N jump to the global maximum/minimum. Hold Shift with
M/N to search within the selected sensor. L switches linear/log frequency
for positive grids and C toggles crosshairs. Interactive navigation requires
a toolkit that delivers pointer and keyboard events (for example, Qt).
With a headless toolkit the figure can still be saved using ``print``.

The Octave figure uses the same map layout and p < .05 mask as the MATLAB
viewer. Phase is pooled on the circle. Color limits are shared across
sensors and observations. The Octave figure exposes standard axes in place
of the App Designer property controls.

References
^^^^^^^^^^
Cohen, M. X. (2014). Analyzing neural time series data: theory and practice. MIT press.
