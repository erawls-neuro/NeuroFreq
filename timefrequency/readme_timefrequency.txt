timefrequency/ README
——————————----------------------------------------------------------------
The timefrequency folder of NeuroFreq contains most of the commonly used TF 
decomposition methods in EEG/neurophysiology research. These functions
include several linear TF decompositions, and a few quadratic TF
distributions. 

Several functions use MATLAB built-in functions for TF transformations. 
These include:
1) nf_stft (stft.m, spectrogram.m, depending on matlab version)
2) nf_cwt (cwt.m)

Several functions use code from GPL software releases, specifically the 
Discrete Time-Frequency Toolbox released by Jeff O'Neill and the 
RID-Rihaczek function released by Selin Aviyente. These functions include:
1) nf_ridbornjordan (O'Neill)
2) nf_ridbinomial (O'Neill)
3) nf_ridrihaczek (Aviyente)

Other functions use custom implementations of TF decomposition methods. 
These include:
1) nf_filterhilbert
2) nf_demodulation
3) nf_dcwt
4) nf_stransform

Many functions use matlab built-in functions other than direct TF 
calculations.

All functions accept input of dimensions 
sensorsXtime OR sensorsXtimeXtrials.

——————————----------------------------------------------------------------
Linear Decompositions: This toolbox provides several linear TF methods.
These methods decompose neural data into frequency band-specific power 
using linear operations.

1) nf_stft:
Short-time Fourier transform using matlab built-in functions. If version
is later than R2019a, uses 'stft' for calculation. If version is earlier
than R2019a, uses 'spectrogram' for calculation. Returns power and phase 
estimates. 

2) nf_filterhilbert: 
Uses a combination of matlab 'butterworth' and 'envelope' functions to 
filter data into specified frequency bands and return the Hilbert envelope 
at those bands. Returns power and phase estimates.

3) nf_demodulation:
Complex demodulation. Multiplies the data with a series of complex 
oscillations (demodulation), then low-pass filters the data to recover the 
demodulated frequency. Returns power and phase estimates.

4) nf_dcwt:
Custom Morlet wavelets described in Cohen (2014). Allows adaptive changes 
in Wavelet time-frequency spread by changing cycle number. Wavelets are 
power-normalized in the frequency domain by dividing the wavelet FFT power 
by its max. Returns power and phase estimates.

5) nf_cwt:
Uses matlab built-in function 'cwt' to compute the continuous wavelet 
transform. CWT returns log-spaced frequencies up to nyquist. Returns power 
and phase estimates.

6) nf_stransform:
Uses the Stockwell transform (S-transform) to compute time-frequency 
distribution. The S-transform is essentially a mix of Fourier and 
wavelet-based methods, modifying the resolution of an output to adapt over 
frequencies using Gaussian windows in frequency domain. Returns power and 
phase estimates.

——————————----------------------------------------------------------------
Quadratic TF Distributions: This toolbox supports several methods 
that involve reflecting the signal times back on itself, prior to 
solving for an entire TF distribution at once. These methods are called
quadratic time-frequency distributions. They can obtain high-resolution 
results, but can also have high interference. Due to the problem of 
interferences in everyday applications the toolbox only supports the 
so-called reduced interference distributions in this class.

7) nf_ridbinomial:
Cohen’s class binomial reduced interference distribution. 
Provides high-resolution TF distributions by windowing the wigner
distribution. Includes inline code from Jeff O’Neill’s dtfd toolbox. 
Returns power but not phase estimates.

8) nf_ridbornjordan:
Cohen’s class Born-Jordan reduced interference distribution. 
Provides high-resolution TF distributions by windowing the wigner
distribution. Includes inline code from Jeff O’Neill’s dtfd toolbox. 
Returns power but not phase estimates.

9) nf_ridrihaczek:
Cohen’s class reduced interference Rihaczek distribution. 
Provides high-resolution complex TF distributions by windowing the 
Rihaczek complex energy spectrum. Includes inline code from Selin 
Aviyente’s RID-Rihaczek toolbox. Returns power and phase estimates.

