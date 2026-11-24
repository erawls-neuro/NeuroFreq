timefrequency/ README
——————————----------------------------------------------------------------
The timefrequency folder of NeuroFreq contains most of the commonly used TF 
decomposition methods in EEG/neurophysiology research. These functions
include several linear TF decompositions, and a few quadratic TF
distributions. 

nf_cwt uses the MATLAB Wavelet Toolbox cwt function.

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
sensorsXtime OR sensorsXtimeXtrials. Each trial is transformed independently.
For multiple sensors, output arrays are sensors x frequencies x times x trials.
With one sensor, output arrays are frequencies x times x trials.

Phase-bearing outputs record their convention in phaseReference.
For frame-start-carrier (STFT) or epoch-start-carrier (Stockwell), add
2*pi*frequency*phaseReferenceSeconds to the coefficient phase to express it
at the reported sample. The offset is in seconds, regardless of EEG event
time zero. Rihaczek phase is the angle of a bilinear distribution.

——————————----------------------------------------------------------------
Linear Decompositions: This toolbox provides several linear TF methods.
These methods decompose neural data into frequency band-specific power 
using linear operations.

1) nf_stft:
Short-time Fourier transform using symmetric Hamming windows and frame-wise
FFTs. Reports squared peak amplitude, including single-gain DC and Nyquist
bins for real input. Phase is referenced to the start of each frame. Output
times are frame centers and Fs is the sampling rate of the frame sequence.

2) nf_filterhilbert: 
Uses forward-backward Butterworth filtering in second-order sections
and a Hilbert transform to obtain the analytic signal in each band. Returns power and phase estimates.

3) nf_demodulation:
Complex demodulation. Multiplies the data with a series of complex 
oscillations (demodulation), then low-pass filters the data to recover the 
demodulated frequency. Returns power and phase estimates.

4) nf_dcwt:
Custom Morlet wavelets described in Cohen (2014). Allows adaptive changes 
in Wavelet time-frequency spread by changing cycle number. Wavelets are 
power-normalized in the frequency domain by dividing the wavelet FFT
by its maximum magnitude. Returns power and phase estimates.

5) nf_cwt:
Computes the analytic Morlet CWT using MATLAB Wavelet Toolbox or the
NeuroFreq FFT implementation in Octave. Both return positive-rotation
power and instantaneous phase on a log-spaced Hz grid. Octave uses ten
voices per octave and symmetric extension within each trial. TF.cwtBackend
identifies the implementation. Compare TF.freqs before pooling runtimes.

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

