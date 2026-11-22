#include "mex.h"
#include <math.h>

/* Positive-lag local autocorrelation for a dense real or complex EEG trace.
 * Build using mex hnf_lacf2_mex.c (legacy and interleaved APIs supported).
 * The MATLAB fallback in hnf_lacf2.m uses the identical lag convention. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1 || nrhs > 2 || nlhs > 1)
        mexErrMsgIdAndTxt("hnf_lacf2_mex:usage", "Use L = hnf_lacf2_mex(x, mlag).");

    const mxArray *x = prhs[0];
    if (!mxIsDouble(x) || mxIsSparse(x) || mxIsEmpty(x) ||
        mxGetNumberOfDimensions(x) != 2 || (mxGetM(x) != 1 && mxGetN(x) != 1))
        mexErrMsgIdAndTxt("hnf_lacf2_mex:x", "x must be a nonempty dense double vector.");

    const mwSize n = mxGetNumberOfElements(x);
    mwSize lags = n;
    if (nrhs == 2 && !mxIsEmpty(prhs[1])) {
        const mxArray *arg = prhs[1];
        if (!mxIsDouble(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
            mexErrMsgIdAndTxt("hnf_lacf2_mex:mlag", "Use a real scalar mlag.");
        const double count = mxGetScalar(arg);
        if (!mxIsFinite(count) || count < 1 || count > (double)n || count != floor(count))
            mexErrMsgIdAndTxt("hnf_lacf2_mex:mlag", "mlag must be an integer in [1, length(x)].");
        lags = (mwSize)count;
    }
    /* Avoid overflow of the element count before requesting the matrix. */
    if ((size_t)lags > ((size_t)-1) / sizeof(double) / (size_t)n / (mxIsComplex(x) ? 2 : 1))
        mexErrMsgIdAndTxt("hnf_lacf2_mex:size", "The requested output exceeds the addressable array size.");
    mxArray *out = mxCreateDoubleMatrix(lags, n, mxIsComplex(x) ? mxCOMPLEX : mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    const mxComplexDouble *z = mxIsComplex(x) ? mxGetComplexDoubles(x) : NULL;
    mxComplexDouble *result = mxIsComplex(x) ? mxGetComplexDoubles(out) : NULL;
    const double *realX = mxIsComplex(x) ? NULL : mxGetDoubles(x);
    double *realOut = mxIsComplex(x) ? NULL : mxGetDoubles(out);
#else
    const double *realX = mxGetPr(x), *imagX = mxGetPi(x);
    double *realOut = mxGetPr(out), *imagOut = mxGetPi(out);
#endif

    for (mwSize t = 0; t < n; ++t) {
        const mwSize count = lags < n-t ? lags : n-t;
        for (mwSize lag = 0; lag < count; ++lag) {
            const mwSize next = t+lag, i = t*lags+lag;
            if (mxIsComplex(x)) {
#if MX_HAS_INTERLEAVED_COMPLEX
                result[i].real = z[t].real*z[next].real + z[t].imag*z[next].imag;
                result[i].imag = z[t].real*z[next].imag - z[t].imag*z[next].real;
#else
                realOut[i] = realX[t]*realX[next] + imagX[t]*imagX[next];
                imagOut[i] = realX[t]*imagX[next] - imagX[t]*realX[next];
#endif
            } else {
                realOut[i] = realX[t]*realX[next];
            }
        }
    }
    if (nlhs) plhs[0] = out;
    else mxDestroyArray(out);
}
