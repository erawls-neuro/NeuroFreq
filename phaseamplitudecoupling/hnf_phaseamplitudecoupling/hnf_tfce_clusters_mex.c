#include "mex.h"
#include <math.h>
#include <float.h>
#include <stdint.h>

/* Integrate signed cluster extent over the same 50 threshold midpoints as
 * hnf_tfce_integrate. The queue holds each active voxel at most once per
 * threshold, even when the montage lists duplicate neighbours. Coordinate
 * and queue indices use mwSize; there are no channel/sample short limits. */
static int real_dense_double(const mxArray *value)
{
    return mxIsDouble(value) && !mxIsComplex(value) && !mxIsSparse(value);
}

void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[])
{
    if (nin != 3 || nout != 1 || !real_dense_double(in[0]) ||
        !real_dense_double(in[1]) || !real_dense_double(in[2]))
        mexErrMsgIdAndTxt("NeuroFreq:TFCEInputs",
                          "Supply a dense real double map, neighbour matrix and two exponents.");
    const mwSize ndim = mxGetNumberOfDimensions(in[0]);
    const mwSize *dims = mxGetDimensions(in[0]);
    const mwSize nc = dims[0], nf = dims[1], nt = ndim == 3 ? dims[2] : 1;
    const mwSize n = mxGetNumberOfElements(in[0]), width = mxGetN(in[1]);
    if (ndim > 3 || mxGetNumberOfDimensions(in[1]) != 2 ||
        mxGetM(in[1]) != nc || mxGetNumberOfElements(in[2]) != 2)
        mexErrMsgIdAndTxt("NeuroFreq:TFCEInputs",
                          "Expected channel-by-sample or channel-by-frequency-by-time with one neighbour row per channel.");
    const double *x = mxGetPr(in[0]), *neighbours = mxGetPr(in[1]), *eh = mxGetPr(in[2]);
    if (!mxIsFinite(eh[0]) || !mxIsFinite(eh[1]) || eh[0] < 0 || eh[1] < 0)
        mexErrMsgIdAndTxt("NeuroFreq:InvalidTFCEExponents", "TFCE exponents must be finite and nonnegative.");
    for (mwSize j = 0; j < mxGetNumberOfElements(in[1]); ++j) {
        double c = neighbours[j];
        if (!mxIsFinite(c) || c < 0 || c > (double)nc || c != floor(c))
            mexErrMsgIdAndTxt("NeuroFreq:NeighbourMismatch", "Invalid one-based channel neighbour.");
    }
    double peak = 0;
    for (mwSize j = 0; j < n; ++j) {
        if (!mxIsFinite(x[j]))
            mexErrMsgIdAndTxt("NeuroFreq:UndefinedTFCE", "TFCE requires finite real statistics.");
        if (fabs(x[j]) > peak) peak = fabs(x[j]);
    }
    out[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
    double *y = mxGetPr(out[0]);
    if (peak == 0) return;
    if ((uintmax_t)n > SIZE_MAX / sizeof(mwSize))
        mexErrMsgIdAndTxt("NeuroFreq:TFCESize", "The voxel queue exceeds this platform's address space.");
    mwSize *queue = (mwSize *)mxMalloc((size_t)n * sizeof(mwSize));
    signed char *active = (signed char *)mxMalloc((size_t)n * sizeof(signed char));
    const double log_peak = log(peak);
    const double extent_bound = eh[0] * log((double)n);
    const double height_bound = eh[1] * fmax(fabs(log_peak), fabs(log_peak - log(100.0)));
    const int use_log = extent_bound + height_bound + fmax(fabs(log_peak),fabs(log_peak-log(50.0))) > -log(DBL_MIN);
    const double delta = (use_log ? 1.0 : peak) / 50.0;
    if (use_log) {
        for (mwSize j = 0; j < n; ++j) y[j] = -mxGetInf();
    }
    for (int step = 0; step < 50; ++step) {
        const double height = delta / 2 + step * delta;
        const double height_weight = use_log ?
            eh[1]*(log_peak + log(height)) + (log_peak-log(50.0)) : pow(height, eh[1])*delta;
        for (mwSize j = 0; j < n; ++j) {
            const double value = use_log ? x[j] / peak : x[j];
            active[j] = value >= height ? 1 : value <= -height ? -1 : 0;
        }
        for (mwSize seed = 0; seed < n; ++seed) {
            const signed char sign = active[seed];
            if (!sign) continue;
            mwSize count = 1;
            active[seed] = 0;
            queue[0] = seed;
            for (mwSize head = 0; head < count; ++head) {
                const mwSize point = queue[head], c = point % nc, ft = point / nc;
                const mwSize f = ft % nf, t = ft / nf;
                const mwSize last_t = t < nt - 1 ? t + 1 : t;
                const mwSize last_f = f < nf - 1 ? f + 1 : f;
                for (mwSize ti = t > 0 ? t - 1 : 0; ti <= last_t; ++ti) {
                    for (mwSize fi = f > 0 ? f - 1 : 0; fi <= last_f; ++fi) {
                        const mwSize offset = nc * (fi + nf * ti);
                        for (mwSize k = 0; k < width; ++k) {
                            const double label = neighbours[c + nc*k];
                            if (label == 0) break;
                            const mwSize near = offset + (mwSize)label - 1;
                            if (active[near] != sign) continue;
                            active[near] = 0;
                            queue[count++] = near;
                        }
                    }
                }
            }
            if (use_log) {
                const double increment = eh[0]*log((double)count) + height_weight;
                for (mwSize j = 0; j < count; ++j) {
                    const mwSize node = queue[j];
                    const double high = fmax(y[node], increment);
                    const double low = fmin(y[node], increment);
                    y[node] = mxIsFinite(high) ? high + log1p(exp(low-high)) : high;
                }
            } else {
                const double increment = sign * (pow((double)count, eh[0]) * height_weight);
                for (mwSize j = 0; j < count; ++j) y[queue[j]] += increment;
            }
        }
    }
    if (use_log) {
        for (mwSize j = 0; j < n; ++j) y[j] = copysign(exp(y[j]), x[j]);
    }
    mxFree(active);
    mxFree(queue);
}
