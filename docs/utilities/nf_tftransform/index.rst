
nf_tftransform
==============

One main feature of the NeuroFreq Toolbox is the inclusion of a high-level function for TF computation, which provides access to every TF algorithm using a unified syntax and a single function with the option to specify the method of choice. The function accepts EEGLAB-formatted (.set) data structures as input. It is called as follows:

.. code-block:: matlab
   
  TF = nf_tftransform(EEG, 'method', 'dcwt', 'freqs', 2:30, 'cycles', [3 8]);

Where 1) EEG is an EEGLAB-formatted (.set) data structure, and 2) the method can be any of the transforms detailed in this package. Additional keyword-argument pairs differ by algorithm.