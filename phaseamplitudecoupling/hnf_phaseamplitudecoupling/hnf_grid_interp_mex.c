/* Linear interpolation of dense real TF traces on a shared increasing axis.
 * Compute the bracket and weight once per target coordinate, reuse them for
 * all columns, and never expand either to a query-by-traces matrix. MATLAB
 * handles phase unwrapping/rewrapping and all nonfinite input data. */
#include "mex.h"

static int realDouble(const mxArray *a)
{
    return mxIsDouble(a) && !mxIsSparse(a) && !mxIsComplex(a);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const double *x, *q, *data;
    double *out, *weight;
    mwSize n, nq, cols, k, col, left = 0;
    mwIndex *bracket;
    if (nrhs != 3 || nlhs != 1 || !realDouble(prhs[0]) ||
        !realDouble(prhs[1]) || !realDouble(prhs[2]) ||
        mxGetNumberOfDimensions(prhs[2]) != 2)
        mexErrMsgIdAndTxt("NeuroFreq:GridArguments", "Expected a real double grid, query and matrix of traces.");
    n = mxGetNumberOfElements(prhs[0]);
    nq = mxGetNumberOfElements(prhs[1]);
    cols = mxGetN(prhs[2]);
    if (n < 2 || mxGetM(prhs[2]) != n)
        mexErrMsgIdAndTxt("NeuroFreq:GridShape", "Each column must match the input grid.");
    x = mxGetPr(prhs[0]); q = mxGetPr(prhs[1]); data = mxGetPr(prhs[2]);
    for (k=0; k<n; ++k)
        if (!mxIsFinite(x[k]) || (k && x[k]<=x[k-1]))
            mexErrMsgIdAndTxt("NeuroFreq:GridAxis", "Input axis must be finite and strictly increasing.");
    bracket = mxCalloc(nq, sizeof(mwIndex));
    weight = mxCalloc(nq, sizeof(double));
    for (k=0; k<nq; ++k) {
        if (!mxIsFinite(q[k]) || (k && q[k] < q[k-1]))
            mexErrMsgIdAndTxt("NeuroFreq:GridQuery", "Queries must be ordered and finite.");
        if (q[k]<x[0] || q[k]>x[n-1]) { bracket[k]=n; continue; }
        while (left+1 < n && x[left+1] <= q[k]) ++left;
        bracket[k]=left;
        if (q[k]!=x[left]) weight[k]=(q[k]-x[left])/(x[left+1]-x[left]);
    }
    plhs[0] = mxCreateDoubleMatrix(nq, cols, mxREAL);
    out = mxGetPr(plhs[0]);
    for (col=0; col<cols; ++col) {
        const double *trace = data+col*n;
        for (k=0; k<nq; ++k) {
            mwIndex i = bracket[k];
            double y;
            if (i==n) y=mxGetNaN();
            else if (weight[k]==0) y=trace[i];
            else y=trace[i]+weight[k]*(trace[i+1]-trace[i]);
            out[col*nq+k]=y;
        }
    }
    mxFree(weight); mxFree(bracket);
}
