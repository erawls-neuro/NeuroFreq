/* Robust division of matching full real single/double arrays (or one scalar).
 * Split out the two precisions so single division rounds exactly as A./B.
 * Invalid operands produce zero. A finite division that overflows may produce
 * Inf, just as in hnf_robdiv.m. MATLAB handles general implicit expansion.
 */
#include "mex.h"
#include <math.h>

static int sameSize(const mxArray *a, const mxArray *b)
{
    mwSize dims = mxGetNumberOfDimensions(a), k;
    const mwSize *sizeA = mxGetDimensions(a), *sizeB = mxGetDimensions(b);
    if (dims != mxGetNumberOfDimensions(b)) return 0;
    for (k = 0; k < dims; ++k)
        if (sizeA[k] != sizeB[k]) return 0;
    return 1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *a, *b, *shape;
    mwSize k, n, stepA, stepB;
    mxClassID precision;

    if (nrhs != 2 || nlhs != 1)
        mexErrMsgIdAndTxt("NeuroFreq:robdiv:arity","Expected two inputs and one output.");
    a = prhs[0];
    b = prhs[1];
    precision = mxGetClassID(a);
    if ((precision != mxDOUBLE_CLASS && precision != mxSINGLE_CLASS) ||
        mxGetClassID(b) != precision || mxIsComplex(a) || mxIsComplex(b) ||
        mxIsSparse(a) || mxIsSparse(b))
        mexErrMsgIdAndTxt("NeuroFreq:robdiv:type","Expected matching full real single or double inputs.");
    stepA = mxGetNumberOfElements(a) != 1;
    stepB = mxGetNumberOfElements(b) != 1;
    if (stepA && stepB && !sameSize(a,b))
        mexErrMsgIdAndTxt("NeuroFreq:robdiv:size","Expected equal shapes or a scalar operand.");
    shape = stepA ? a : b;
    n = mxGetNumberOfElements(shape);
    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(shape), mxGetDimensions(shape), precision, mxREAL);

    if (precision == mxDOUBLE_CLASS) {
        const double *x = mxGetPr(a), *y = mxGetPr(b);
        double *out = mxGetPr(plhs[0]);
        for (k = 0; k < n; ++k) {
            double u = x[k*stepA], v = y[k*stepB];
            out[k] = (v != 0.0 && isfinite(u) && isfinite(v)) ? u / v : 0.0;
        }
    } else {
        const float *x = (const float *)mxGetData(a), *y = (const float *)mxGetData(b);
        float *out = (float *)mxGetData(plhs[0]);
        for (k = 0; k < n; ++k) {
            float u = x[k*stepA], v = y[k*stepB];
            out[k] = (v != 0.0f && isfinite(u) && isfinite(v)) ? u / v : 0.0f;
        }
    }
}
