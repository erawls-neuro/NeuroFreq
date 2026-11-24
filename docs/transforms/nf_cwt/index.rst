
nf_cwt
======

The Continuous Wavelet Transform (CWT) transforms data to a TF representation using wavelets, which are complex sine waves multiplied by a window function that localizes the wavelet energy at the center time. The CWT scales and translates a “mother wavelet” by dilating it in time, localizing different frequencies in the process. CWT is calculated using

.. code-block:: matlab
   
  TF = nf_cwt( data, Fs, verbose );

Where data is a 1/2/3D tensor of dimensions channels X time X trials, Fs is the sampling rate of the data in Hz, and verbose is 0 or 1 indicating whether to print progress during transformation. Note that for Morlet wavelets, the MATLAB function cwt.m does not allow modifying the number of cycles.

Defaults are: verbose = 0. Data and Fs are required.

When CWT is run on the demo synthetic data using

.. code-block:: matlab
  
  TF = nf_cwt( data, 500 );

We obtain the following result:

.. image:: fig_cwt_synthetic.png
  :width: 600

GNU Octave uses NeuroFreq's FFT-based analytic Morlet bank. Each trial is
reflected at its endpoints before convolution. The bank has ten voices per
octave and a 10% response at Nyquist in its highest band. Frequencies are
returned in ascending Hz. For a complex signal the returned phase follows
positive rotation. The output field ``cwtBackend`` identifies the
implementation. Before combining MATLAB and Octave runs, compare
``TF.freqs`` and check agreement on a representative epoch, particularly
near the edges. The Octave implementation uses an analytic estimate of
wavelet spread to choose its lowest frequency.
