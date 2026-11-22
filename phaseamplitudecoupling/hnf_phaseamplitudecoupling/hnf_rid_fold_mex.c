#include "mex.h"
#include <math.h>
#include <stdint.h>

/* Multiply ambiguity by the real Choi-Williams kernel and accumulate its
 * signed lag rows modulo nfft. Walk contiguous columns in both inputs.
 * Lag ordering and the order of additions match nf_ridrihaczek.m. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3 || nlhs != 1)
        mexErrMsgIdAndTxt("NeuroFreq:RIDFoldArgs","Use F = hnf_rid_fold_mex(ambiguity,kernel,nfft).");
    const mxArray *a = prhs[0], *kernel = prhs[1], *length = prhs[2];
    size_t n = mxGetM(a);
    if (!mxIsDouble(a) || mxIsSparse(a) || n < 1 || mxGetNumberOfDimensions(a)!=2 ||
        mxGetN(a)!=n || !mxIsDouble(kernel) || mxIsSparse(kernel) || mxIsComplex(kernel) ||
        mxGetNumberOfDimensions(kernel)!=2 || mxGetM(kernel)!=n || mxGetN(kernel)!=n ||
        !mxIsDouble(length) || mxIsSparse(length) || mxIsComplex(length) ||
        mxGetNumberOfElements(length)!=1)
        mexErrMsgIdAndTxt("NeuroFreq:RIDFoldInput","Expected two full double square matrices of equal size and a transform length.");
    double count = mxGetScalar(length);
    if (!mxIsFinite(count) || count < 1 || floor(count)!=count ||
        count >= (double)(SIZE_MAX/sizeof(double)/2/(size_t)n))
        mexErrMsgIdAndTxt("NeuroFreq:RIDFoldSize","Invalid FFT size.");
    size_t m = (size_t)count;
    size_t *rows = mxCalloc(n,sizeof(*rows));
    for (size_t lag = 0; lag < n; ++lag) {
        /* bins = ifftshift(-floor(n/2):ceil(n/2)-1) */
        if (lag <= (n-1)/2) rows[lag] = (size_t)lag % m;
        else rows[lag] = (m-((size_t)(n-lag) % m)) % m;
    }
    const double *weights = mxGetPr(kernel);
    int complexA = mxIsComplex(a);
    plhs[0] = mxCreateDoubleMatrix(m,n,complexA ? mxCOMPLEX : mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
    const mxComplexDouble *az = complexA ? mxGetComplexDoubles(a) : NULL;
    mxComplexDouble *fz = complexA ? mxGetComplexDoubles(plhs[0]) : NULL;
    const double *ar = complexA ? NULL : mxGetDoubles(a);
    double *fr = complexA ? NULL : mxGetDoubles(plhs[0]);
#else
    const double *ar = mxGetPr(a), *ai = mxGetPi(a);
    double *fr = mxGetPr(plhs[0]), *fi = mxGetPi(plhs[0]);
#endif
    for (size_t col = 0; col < n; ++col) {
        for (size_t lag = 0; lag < n; ++lag) {
            size_t i = (size_t)lag+(size_t)col*n, dst = rows[lag]+(size_t)col*m;
#if MX_HAS_INTERLEAVED_COMPLEX
            if (complexA) {
                fz[dst].real += az[i].real*weights[i];
                fz[dst].imag += az[i].imag*weights[i];
            } else fr[dst] += ar[i]*weights[i];
#else
            fr[dst] += ar[i]*weights[i];
            if (complexA) fi[dst] += ai[i]*weights[i];
#endif
        }
    }
    mxFree(rows);
}
