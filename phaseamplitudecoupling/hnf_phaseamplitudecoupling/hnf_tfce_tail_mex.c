#include "mex.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Spatially corrected Monte Carlo p-values for one effect.
 * The null vector is copied and sorted once. A binary search on this vector
 * counts exceedances for each voxel without sorting a large observed map.
 * Follow hnf_tfce_pvalues for NaN, signed observations and roundoff ties.
 */
static int compare_double(const void *first, const void *second)
{
    double a = *(const double *)first;
    double b = *(const double *)second;
    return (a > b) - (a < b);
}

static int dense_real_double(const mxArray *x)
{
    return mxIsDouble(x) && !mxIsSparse(x) && !mxIsComplex(x);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nlhs != 1 || nrhs != 2 || !dense_real_double(prhs[0]) || !dense_real_double(prhs[1])) {
        mexErrMsgIdAndTxt("NeuroFreq:TailInputs", "Use p = hnf_tfce_tail_mex(observed, maxima) with full real doubles.");
    }
    const double *observed = mxGetPr(prhs[0]);
    const double *maxima = mxGetPr(prhs[1]);
    mwSize m = mxGetNumberOfElements(prhs[0]);
    mwSize n = mxGetNumberOfElements(prhs[1]);
    if (n == 0) mexErrMsgIdAndTxt("NeuroFreq:InvalidNull", "An empty permutation distribution is undefined.");
    for (mwIndex i = 0; i < n; ++i) {
        if (!mxIsFinite(maxima[i])) mexErrMsgIdAndTxt("NeuroFreq:InvalidNull", "Permutation maxima must be finite.");
    }
    double *sorted = (double *)mxMalloc(n * sizeof(double));
    memcpy(sorted, maxima, n * sizeof(double));
    qsort(sorted, n, sizeof(double), compare_double);
    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]), mxGetDimensions(prhs[0]), mxDOUBLE_CLASS, mxREAL);
    double *p = mxGetPr(plhs[0]);
    for (mwIndex i = 0; i < m; ++i) {
        if (!mxIsFinite(observed[i])) {
            p[i] = mxGetNaN();
            continue;
        }
        double cutoff = fabs(observed[i]) * (1.0 - 100.0 * DBL_EPSILON);
        mwIndex low = 0, high = n;
        while (low < high) {
            mwIndex mid = low + (high - low) / 2;
            if (sorted[mid] < cutoff) low = mid + 1;
            else high = mid;
        }
        p[i] = (1.0 + (double)(n - low)) / (1.0 + (double)n);
    }
    mxFree(sorted);
}
