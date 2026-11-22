#include "mex.h"
#include <math.h>
#include <stdint.h>

/* Pack tapered windows as contiguous columns and run a batch of MATLAB's
 * ordinary FFTs. No additional FFT library is required. Frames are bounded
 * in batches to keep large, highly-overlapping EEG analyses memory-safe.
 * x is [sensor, sample, trial]; outputs are [sensor, bin, frame, trial]. */
static size_t count_arg(const mxArray *a)
{
    if (!mxIsDouble(a) || mxIsSparse(a) || mxIsComplex(a) ||
        mxGetNumberOfElements(a) != 1)
        mexErrMsgIdAndTxt("NeuroFreq:STFTSize","Expected an integer sample count.");
    double x = mxGetScalar(a);
    if (!(x >= 1 && x == floor(x) && x < (double)(SIZE_MAX/sizeof(double)/2)))
        mexErrMsgIdAndTxt("NeuroFreq:STFTSize","Invalid sample count.");
    return (size_t)x;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs!=5 || nlhs!=2)
        mexErrMsgIdAndTxt("NeuroFreq:STFTArgs","Use [power,phase]=hnf_stft_mex(x,taper,nfft,hop,gain).");
    const mxArray *x = prhs[0], *w = prhs[1], *g = prhs[4];
    if (!mxIsDouble(x) || mxIsEmpty(x) || mxIsSparse(x) || mxGetNumberOfDimensions(x)>3 ||
        !mxIsDouble(w) || mxIsSparse(w) || mxIsComplex(w) ||
        !mxIsDouble(g) || mxIsSparse(g) || mxIsComplex(g))
        mexErrMsgIdAndTxt("NeuroFreq:STFTInput","Expected full double data and real window/bin gains.");
    const mwSize *dims = mxGetDimensions(x);
    size_t sensors = dims[0], times = dims[1], trials = mxGetNumberOfDimensions(x)==3 ? (size_t)dims[2] : 1;
    size_t nfft = count_arg(prhs[2]), hop = count_arg(prhs[3]), window = mxGetNumberOfElements(w);
    size_t bins = nfft/2+1;
    if (window<2 || window>times || window>nfft || mxGetNumberOfElements(g)!=bins)
        mexErrMsgIdAndTxt("NeuroFreq:STFTGrid","FFT, taper, and gain sizes do not match the data.");
    size_t frames = 1 + (times-window)/hop;
    if (frames > SIZE_MAX/sensors/trials || bins > SIZE_MAX/sizeof(double)/sensors/trials/frames)
        mexErrMsgIdAndTxt("NeuroFreq:STFTSize","The requested spectrogram is too large.");
    size_t columns = sensors*trials*frames;
    mwSize outDims[4] = {(mwSize)sensors,(mwSize)bins,(mwSize)frames,(mwSize)trials};
    plhs[0] = mxCreateNumericArray(4,outDims,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(4,outDims,mxDOUBLE_CLASS,mxREAL);
    double *power = mxGetPr(plhs[0]), *phase = mxGetPr(plhs[1]);
    const double *taper = mxGetPr(w), *gain = mxGetPr(g);
    int complexX = mxIsComplex(x);
#if MX_HAS_INTERLEAVED_COMPLEX
    const mxComplexDouble *xz = complexX ? mxGetComplexDoubles(x) : NULL;
    const double *xr = complexX ? NULL : mxGetDoubles(x);
#else
    const double *xr = mxGetPr(x), *xi = mxGetPi(x);
#endif
    /* Allow at most ~4 MiB per complex work tensor, or one long frame. */
    size_t batch = (4*1024*1024)/(2*sizeof(double)*nfft);
    if (!batch) batch = 1;
    for (size_t start=0; start<columns;) {
        size_t count = columns-start < batch ? columns-start : batch;
        mxArray *packed = mxCreateDoubleMatrix(nfft,count,complexX ? mxCOMPLEX : mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble *pz = complexX ? mxGetComplexDoubles(packed) : NULL;
        double *pr = complexX ? NULL : mxGetDoubles(packed);
#else
        double *pr = mxGetPr(packed), *pi = mxGetPi(packed);
#endif
        for (size_t col=0; col<count; ++col) {
            size_t job = start+col, sensor = job%sensors, frame = job/sensors%frames, trial = job/sensors/frames;
            size_t base = sensor+sensors*(frame*hop+trial*times);
            for (size_t sample=0; sample<window; ++sample) {
                size_t dst = sample+nfft*col, src = base+sample*sensors;
#if MX_HAS_INTERLEAVED_COMPLEX
                if (complexX) {
                    pz[dst].real = taper[sample]*xz[src].real;
                    pz[dst].imag = taper[sample]*xz[src].imag;
                } else pr[dst] = taper[sample]*xr[src];
#else
                pr[dst] = taper[sample]*xr[src];
                if (complexX) pi[dst] = taper[sample]*xi[src];
#endif
            }
        }
        mxArray *spectrum = NULL;
        if (mexCallMATLAB(1,&spectrum,1,&packed,"fft")!=0)
            mexErrMsgIdAndTxt("NeuroFreq:STFTFFT","The MATLAB frame FFT failed.");
        mxDestroyArray(packed);
        int complexS = mxIsComplex(spectrum);
#if MX_HAS_INTERLEAVED_COMPLEX
        const mxComplexDouble *sz = complexS ? mxGetComplexDoubles(spectrum) : NULL;
        const double *sr = complexS ? NULL : mxGetDoubles(spectrum);
#else
        const double *sr = mxGetPr(spectrum), *si = mxGetPi(spectrum);
#endif
        for (size_t col=0; col<count; ++col) {
            size_t job = start+col, sensor = job%sensors;
            size_t base = sensor + (job/sensors)*sensors*bins;
            for (size_t k=0; k<bins; ++k) {
                size_t src = col*nfft+k, dst = base+k*sensors;
#if MX_HAS_INTERLEAVED_COMPLEX
                double re = gain[k]*(complexS ? sz[src].real : sr[src]);
                double im = complexS ? gain[k]*sz[src].imag : 0.0;
#else
                double re = gain[k]*sr[src], im = complexS ? gain[k]*si[src] : 0.0;
#endif
                double magnitude = hypot(re,im);
                power[dst] = magnitude*magnitude;
                phase[dst] = atan2(im,re);
            }
        }
        mxDestroyArray(spectrum);
        start += count;
    }
}
