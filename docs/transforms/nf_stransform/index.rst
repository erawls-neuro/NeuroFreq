
nf_stransform
=============

The S-transform or Stockwell transform  combines the advantages of both STFT and wavelet methods, and can potentially provide a clearer, more detailed time-frequency representation of M/EEG data compared to either method used alone. The S-transform is calculated using 

.. code-block:: matlab
   
  TF = nf_stransform( data, Fs, factor, verbose );

Where data is a 1/2/3D tensor of dimensions channels X time X trials, Fs is the sampling rate of the data in Hz, factor scales the frequency-dependent Gaussian window, and verbose is 0 or 1 indicating whether to print progress during transformation.

The Gaussian scale factor defaults to 1 and verbose defaults to 0. Data and Fs are required.

When S-transform is run on the demo synthetic data using

.. code-block:: matlab
  
  TF = nf_stransform( data, 500 );

We obtain the following result:

.. image:: fig_stransform_synthetic.png
  :width: 600

Note that the S-transform has excellent resolution at low-frequencies.

References
^^^^^^^^^^

Stockwell, R. G., Mansinha, L., & Lowe, R. P. (1996). Localization of the complex spectrum: the S transform. IEEE transactions on signal processing, 44(4), 998-1001.