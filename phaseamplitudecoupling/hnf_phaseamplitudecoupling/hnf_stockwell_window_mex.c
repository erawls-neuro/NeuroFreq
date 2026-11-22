#include "mex.h"
#include <math.h>

/* Each column k carries the same frequency-shifted, Gaussian-weighted
 * spectrum used by the scalar Stockwell implementation. The caller runs
 * MATLAB's ifft down all columns together. Real and imaginary sums never
 * mix, so both MEX complex-data APIs support single/double input. */
#define WINDOW_SEPARATE(T,src,dst) \
    for (size_t col=0; col<count; ++col) { \
        for (size_t t=0; t<n; ++t) { \
            size_t i=t+col*n; (dst)[i]=(src)[t+first+col]*(T)g[i]; \
        } \
    }
#define WINDOW_INTERLEAVED(T,src,dst) \
    for (size_t col=0; col<count; ++col) { \
        for (size_t t=0; t<n; ++t) { \
            size_t i=t+col*n, j=t+first+col; T gain=(T)g[i]; \
            (dst)[i].real=(src)[j].real*gain; (dst)[i].imag=(src)[j].imag*gain; \
        } \
    }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs!=3 || nlhs!=1)
        mexErrMsgIdAndTxt("NeuroFreq:STWindowArgs","Use block = hnf_stockwell_window_mex(doubledSpectrum,gaussians,firstBin).");
    const mxArray *x=prhs[0], *windows=prhs[1], *firstArg=prhs[2];
    size_t n=mxGetM(windows), count=mxGetN(windows);
    if (!(mxIsDouble(x)||mxIsSingle(x)) || mxIsSparse(x) || mxGetNumberOfDimensions(x)!=2 ||
        (mxGetN(x)!=1 && mxGetM(x)!=1) || mxGetNumberOfElements(x)/2!=n || mxGetNumberOfElements(x)%2 ||
        !mxIsDouble(windows) || mxIsSparse(windows) || mxIsComplex(windows) ||
        mxGetNumberOfDimensions(windows)!=2 || n<2 || !count || count>n/2 ||
        !mxIsDouble(firstArg) || mxIsSparse(firstArg) || mxIsComplex(firstArg) ||
        mxGetNumberOfElements(firstArg)!=1)
        mexErrMsgIdAndTxt("NeuroFreq:STWindowInput","Expected a repeated FFT vector, N-by-B real Gaussian weights and the starting positive FFT bin.");
    double index=mxGetScalar(firstArg);
    if (!(index>=1 && index==floor(index) && index<=(double)(n/2-count+1)))
        mexErrMsgIdAndTxt("NeuroFreq:STWindowBin","The requested bins must lie within the positive DFT frequencies.");
    size_t first=(size_t)index;
    const double *g=mxGetPr(windows);
    int isComplex=mxIsComplex(x);
    plhs[0]=mxCreateNumericMatrix(n,count,mxGetClassID(x),isComplex ? mxCOMPLEX : mxREAL);
    if (mxIsDouble(x)) {
#if MX_HAS_INTERLEAVED_COMPLEX
        if (isComplex) {
            const mxComplexDouble *src=mxGetComplexDoubles(x);
            mxComplexDouble *dst=mxGetComplexDoubles(plhs[0]);
            WINDOW_INTERLEAVED(double,src,dst);
        } else {
            const double *src=mxGetDoubles(x); double *dst=mxGetDoubles(plhs[0]);
            WINDOW_SEPARATE(double,src,dst);
        }
#else
        const double *src=mxGetPr(x); double *dst=mxGetPr(plhs[0]);
        WINDOW_SEPARATE(double,src,dst);
        if (isComplex) {
            src=mxGetPi(x); dst=mxGetPi(plhs[0]);
            WINDOW_SEPARATE(double,src,dst);
        }
#endif
    } else {
#if MX_HAS_INTERLEAVED_COMPLEX
        if (isComplex) {
            /* Octave 9's typed interleaved-single accessor rejects a valid
             * complex single. Its generic accessor uses the same ABI layout. */
#ifdef HAVE_OCTAVE
            const mxComplexSingle *src=(const mxComplexSingle *)mxGetData(x);
            mxComplexSingle *dst=(mxComplexSingle *)mxGetData(plhs[0]);
#else
            const mxComplexSingle *src=mxGetComplexSingles(x);
            mxComplexSingle *dst=mxGetComplexSingles(plhs[0]);
#endif
            WINDOW_INTERLEAVED(float,src,dst);
        } else {
            const float *src=mxGetSingles(x); float *dst=mxGetSingles(plhs[0]);
            WINDOW_SEPARATE(float,src,dst);
        }
#else
        const float *src=mxGetData(x); float *dst=mxGetData(plhs[0]);
        WINDOW_SEPARATE(float,src,dst);
        if (isComplex) {
            src=mxGetImagData(x); dst=mxGetImagData(plhs[0]);
            WINDOW_SEPARATE(float,src,dst);
        }
#endif
    }
}
