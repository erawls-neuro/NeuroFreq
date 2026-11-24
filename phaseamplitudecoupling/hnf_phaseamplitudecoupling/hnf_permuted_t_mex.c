#include "mex.h"
#include <math.h>

/* One observation-order-preserving two-pass t statistic per spatial voxel.
 * The caller generates the same signs or group assignments as the MATLAB
 * path. Assemble no observations-by-voxels permuted array and draw no RNG.
 *
 * hnf_permuted_t_mex(D, signs)           sign flip, single-sample t
 * hnf_permuted_t_mex(D, rowOrder, nA)   group reassignment, Welch t
 *
 * First dimension is observations. Return [1, spatial dimensions...].
 */
static int dense_real_double(const mxArray *x)
{
    return mxIsDouble(x) && !mxIsComplex(x) && !mxIsSparse(x);
}

/* Determine a common scale for a voxel before squaring. Keep ordinary
 * units on the fast path. Powers of two rescale without rounding the
 * significand, including inputs near DBL_MAX or the subnormal range. */
static double safe_scale(const double *x, mwSize n)
{
    double maximum = 0.0;
    for (mwIndex i = 0; i < n; ++i) {
        double value = fabs(x[i]);
        if (value > maximum) maximum = value;
    }
    if (maximum == 0.0 || (maximum < 1e120 && maximum > 1e-120)) return 1.0;
    int exponent;
    frexp(maximum, &exponent);
    return ldexp(1.0, exponent - 1);
}

static double one_sample(const double *x, const double *signs, mwSize n)
{
    double origin = signs[0] * x[0];
    double meanDelta = 0.0, squares = 0.0, total = 0.0, compensation = 0.0;
    for (mwIndex i = 0; i < n; ++i) {
        double term = signs[i] * x[i], next = total + term;
        compensation += fabs(total) >= fabs(term) ? (total-next)+term : (term-next)+total;
        total = next;
        meanDelta += term - origin;
    }
    meanDelta /= (double)n;
    for (mwIndex i = 0; i < n; ++i) {
        double delta = (signs[i] * x[i] - origin) - meanDelta;
        squares += delta * delta;
    }
    return ((total+compensation)/(double)n) / (sqrt(squares / (double)(n-1)) / sqrt((double)n));
}

static double group_mean_delta(const double *x, const mwIndex *order,
                               mwSize n, double origin, double *se)
{
    double mean = 0.0, squares = 0.0, maximum = 0.0, scale = 1.0;
    for (mwIndex i = 0; i < n; ++i) {
        double delta = x[order[i]] - origin;
        mean += delta;
        if (fabs(delta) > maximum) maximum = fabs(delta);
    }
    /* The group with small dispersion can carry the entire standard error.
     * Rescale its residuals separately before squaring them. */
    if (maximum > 0.0 && maximum < 1e-120) {
        int exponent;
        frexp(maximum, &exponent);
        scale = ldexp(1.0, exponent - 1);
        mean = 0.0;
        for (mwIndex i = 0; i < n; ++i) mean += (x[order[i]]-origin)/scale;
    }
    mean /= (double)n;
    for (mwIndex i = 0; i < n; ++i) {
        double centered = (x[order[i]] - origin)/scale - mean;
        squares += centered * centered;
    }
    *se = (sqrt(squares / (double)(n - 1)) / sqrt((double)n))*scale;
    return mean*scale;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nlhs != 1 || (nrhs != 2 && nrhs != 3) || !dense_real_double(prhs[0]) || !dense_real_double(prhs[1])) {
        mexErrMsgIdAndTxt("NeuroFreq:PermutationInputs", "Use t = hnf_permuted_t_mex(data, signs) or t = hnf_permuted_t_mex(data, rowOrder, nA) with full real doubles.");
    }
    mwSize n = mxGetM(prhs[0]);
    if (n < 2 || mxGetNumberOfElements(prhs[1]) != n || (mxGetM(prhs[1]) != 1 && mxGetN(prhs[1]) != 1)) {
        mexErrMsgIdAndTxt("NeuroFreq:PermutationSize", "Supply one sign or one row index for every observation.");
    }
    const double *data = mxGetPr(prhs[0]);
    const double *draw = mxGetPr(prhs[1]);
    mwSize nA = 0;
    mwIndex *order = NULL;
    if (nrhs == 2) {
        for (mwIndex i = 0; i < n; ++i) {
            if (draw[i] != -1.0 && draw[i] != 1.0) mexErrMsgIdAndTxt("NeuroFreq:InvalidSigns", "Signs must be -1 or +1.");
        }
    } else {
        if (!dense_real_double(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1) {
            mexErrMsgIdAndTxt("NeuroFreq:InvalidGroup", "nA must be a scalar integer.");
        }
        double group = mxGetScalar(prhs[2]);
        if (!mxIsFinite(group) || group != floor(group) || group < 2 || group > (double)n - 2) {
            mexErrMsgIdAndTxt("NeuroFreq:InvalidGroup", "Each group must contain at least two observations.");
        }
        nA = (mwSize)group;
        order = (mwIndex *)mxMalloc(n * sizeof(mwIndex));
        mxLogical *seen = (mxLogical *)mxCalloc(n, sizeof(mxLogical));
        for (mwIndex i = 0; i < n; ++i) {
            double value = draw[i];
            if (!mxIsFinite(value) || value < 1 || value > (double)n || value != floor(value)) {
                mexErrMsgIdAndTxt("NeuroFreq:InvalidIndex", "rowOrder must contain valid integer row indices.");
            }
            order[i] = (mwIndex)value - 1;
            if (seen[order[i]]) mexErrMsgIdAndTxt("NeuroFreq:RepeatedIndex", "rowOrder must be a permutation.");
            seen[order[i]] = 1;
        }
        mxFree(seen);
    }

    mwSize dimensions = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *inputSize = mxGetDimensions(prhs[0]);
    mwSize *outputSize = (mwSize *)mxMalloc(dimensions * sizeof(mwSize));
    for (mwIndex k = 0; k < dimensions; ++k) outputSize[k] = inputSize[k];
    outputSize[0] = 1;
    plhs[0] = mxCreateNumericArray(dimensions, outputSize, mxDOUBLE_CLASS, mxREAL);
    mxFree(outputSize);
    double *result = mxGetPr(plhs[0]);
    mwSize voxels = mxGetNumberOfElements(plhs[0]);
    double *scaled = NULL;
    for (mwIndex j = 0; j < voxels; ++j) {
        const double *x = data + j*n;
        double scale = safe_scale(x, n);
        if (scale != 1.0) {
            if (!scaled) scaled = (double *)mxMalloc(n * sizeof(double));
            for (mwIndex i = 0; i < n; ++i) scaled[i] = x[i] / scale;
            x = scaled;
        }
        if (nrhs == 2) {
            result[j] = one_sample(x, draw, n);
        } else {
            double seA, seB, origA = x[order[0]], origB = x[order[nA]];
            double da = group_mean_delta(x, order, nA, origA, &seA);
            double db = group_mean_delta(x, order+nA, n-nA, origB, &seB);
            result[j] = ((origA-origB)+(da-db)) / hypot(seA,seB);
        }
    }
    if (scaled != NULL) mxFree(scaled);
    if (order != NULL) mxFree(order);
}
