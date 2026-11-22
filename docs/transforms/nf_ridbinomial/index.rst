
nf_ridbinomial
==============

The Binomial Reduced Interference Distribution, a member of Cohen's (1995) class, adopts a cross-shaped kernel in the ambiguity domain (Jeong & Williams, 1992). The binomial RID is computed using

.. code-block:: matlab
   
  TF = nf_ridbinomial( data, Fs, fRes, kernel, makePos, verbose );

where data is a 1/2/3D tensor of dimensions channels X time X trials, Fs is the sampling rate of the data in Hz, fRes is the desired frequency resolution of the output in Hz (e.g. 0.5 produces a TF representation with frequency steps equal to 0.5 Hz), kernel is the length of the binomial kernel in seconds, makePos is 0 or 1 indicating whether the distribution should return only positive energy, and verbose is 0 or 1 indicating whether to print progress during transformation. Shorter kernels more effectively mask cross-terms, but also smear the TF representation. The distribution is signed. A positive display floor changes spectral values and their marginals; leave makePos at zero when analyzing the signed distribution.

Defaults are: fRes = as many frequency samples as time samples, kernel = the largest even lag window permitted by the sample and frequency grids, makePos = 0, verbose = 0. Data and Fs are required.

When the binomial RID is run on the demo synthetic data using

.. code-block:: matlab
  
  TF = nf_ridbinomial( data, 500 );

We obtain the following result:

.. image:: fig_binomial_synthetic.png
  :width: 600

References
^^^^^^^^^^
Cohen, L. (1995). Time-frequency analysis (Vol. 778). New Jersey: Prentice hall.

Jeong, J., & Williams, W. J. (1992). Kernel design for reduced interference distributions. IEEE Transactions on Signal Processing, 40(2), 402-412.