/* Finite real columns chosen for compensation by hnf_robsum. Scaling and
 * Neumaier compensation match elem_compsum. Exceptional exponent spreads
 * and complex-component extraction are handled by the MATLAB caller. */
#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize n, columns, c, k;
    const double *x, *scale;
    double *out;
    if (nrhs != 2 || nlhs != 1 ||
        !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) ||
        !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) ||
        mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgIdAndTxt("NeuroFreq:CompensatedSum", "Supply a real double matrix and one scale per column.");
    n = mxGetM(prhs[0]); columns = mxGetN(prhs[0]);
    if (mxGetNumberOfElements(prhs[1]) != columns)
        mexErrMsgIdAndTxt("NeuroFreq:CompensatedSum", "Each column requires a scale.");
    x = mxGetPr(prhs[0]); scale = mxGetPr(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(1,columns,mxREAL); out = mxGetPr(plhs[0]);
    for (c = 0; c < columns; ++c) {
        double total = 0, correction = 0;
        if (!(scale[c]>0 && mxIsFinite(scale[c])))
            mexErrMsgIdAndTxt("NeuroFreq:CompensatedSum", "Scales must be positive and finite.");
        for (k = 0; k < n; ++k) {
            double a = x[c*n+k]/scale[c], next = total+a;
            correction += fabs(total)>=fabs(a) ? (total-next)+a : (a-next)+total;
            total = next;
        }
        out[c] = (total+correction)*scale[c];
    }
}
