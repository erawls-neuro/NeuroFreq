
nf_demodulation
===============

Complex Demodulation calculates time-frequency representations by first multiplying the real-valued signal by a complex sine wave of a specific frequency, which demodulates the power at that frequency so that it is concentrated at low frequencies. Then, the signal is low-pass filtered in order to recover power and phase at the frequency of the sine wave. Complex demodulation is the TF method of choice in BESA software (Hoechstetter et al., 2004). It is calculated using

.. code-block:: matlab
   
  TF = nf_demodulation( data, Fs, freqs, lowpassF, order, verbose );

Where data is a 1/2/3D tensor of dimensions channels X time X trials, Fs is the sampling rate of the data in Hz, freqs is a vector of center frequencies for the sine waves, lowpassF is the frequency of the low-pass filter in Hz, order is the order of the low-pass filter, and verbose is 0 or 1 indicating whether to print progress during transformation.

Defaults are: freqs = 1:1:Fs/2 Hz, lowpassF = 2 Hz, order = 3, verbose = 0. Data and Fs are required.

When demodulation is run on the demo synthetic data using

.. code-block:: matlab
  
  TF = nf_demodulation( data, 500 );

We obtain the following result:

.. image:: fig_demodulation_synthetic.png
  :width: 600


References
^^^^^^^^^^
Hoechstetter, K., Bornfleth, H., Weckesser, D., Ille, N., Berg, P., & Scherg, M. (2004). BESA source coherence: a new method to study cortical oscillatory coupling. Brain topography, 16, 233-238.
