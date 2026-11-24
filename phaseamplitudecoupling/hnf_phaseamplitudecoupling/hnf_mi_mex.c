#include "mex.h"
#include <math.h>

/* Tort MI over the first axis. Reuse one pair of bin arrays across series.
 * The MATLAB fallback documents the missing-bin and invalid-weight policy. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3 || nlhs > 1)
        mexErrMsgIdAndTxt("NeuroFreq:MIArguments", "Use mi = hnf_mi_mex(phase, weights, nbin).");
    for (int k = 0; k < 3; ++k)
        if (!mxIsDouble(prhs[k]) || mxIsSparse(prhs[k]) || mxIsComplex(prhs[k]))
            mexErrMsgIdAndTxt("NeuroFreq:MIInput", "Inputs must be full, real double arrays.");
    mwSize nd = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *dims = mxGetDimensions(prhs[0]);
    if (mxIsEmpty(prhs[0]) || mxGetNumberOfDimensions(prhs[1]) != nd)
        mexErrMsgIdAndTxt("NeuroFreq:MIInput", "Nonempty arrays of equal size are required.");
    for (mwSize k = 0; k < nd; ++k)
        if (dims[k] != mxGetDimensions(prhs[1])[k])
            mexErrMsgIdAndTxt("NeuroFreq:MIInput", "Phase and weights must have equal size.");
    double nbin = mxGetScalar(prhs[2]);
    if (mxGetNumberOfElements(prhs[2]) != 1 || !mxIsFinite(nbin) || nbin < 2 || floor(nbin) != nbin)
        mexErrMsgIdAndTxt("NeuroFreq:MIBins", "nbin must be a finite integer >= 2.");

    mwSize samples = dims[0], series = mxGetNumberOfElements(prhs[0]) / samples;
    mwSize *shape = mxCalloc(nd, sizeof(*shape));
    for (mwSize k = 0; k < nd; ++k) shape[k] = k ? dims[k] : 1;
    plhs[0] = mxCreateNumericArray(nd, shape, mxDOUBLE_CLASS, mxREAL);
    mxFree(shape);
    const double *phase = mxGetPr(prhs[0]), *weights = mxGetPr(prhs[1]);
    double *result = mxGetPr(plhs[0]);
    for (mwSize v = 0; v < series; ++v) result[v] = mxGetNaN();
    if (nbin > (double)samples) return;
    mwSize bins = (mwSize)nbin;
    double *means = mxCalloc(bins, sizeof(*means));
    mwSize *counts = mxCalloc(bins, sizeof(*counts));
    const double pi = 3.14159265358979323846, period = 2*pi;
    for (mwSize v = 0; v < series; ++v) {
        const double *x = phase + v*samples, *w = weights + v*samples;
        double scale = 0, total = 0, kl = 0;
        for (mwSize t = 0; t < samples; ++t) {
            if (!mxIsFinite(x[t]) || !mxIsFinite(w[t]) || w[t] < 0) { scale = 0; break; }
            if (w[t] > scale) scale = w[t];
        }
        if (scale == 0) continue;
        for (mwSize b = 0; b < bins; ++b) { means[b] = 0; counts[b] = 0; }
        for (mwSize t = 0; t < samples; ++t) {
            double theta = x[t];
            if (fabs(theta)>1048576.0) theta = atan2(sin(theta),cos(theta));
            double angle = fmod(theta+pi, period);
            if (angle < 0) angle += period;
            mwSize b = (mwSize)floor(angle/(period/nbin));
            if (b >= bins) b = bins-1;
            means[b] += w[t]/scale;
            ++counts[b];
        }
        for (mwSize b = 0; b < bins; ++b) {
            if (counts[b] == 0) { total = 0; break; }
            means[b] /= (double)counts[b];
            total += means[b];
        }
        if (total == 0) continue;
        for (mwSize b = 0; b < bins; ++b) {
            double p = means[b]/total;
            if (p > 0) kl += p*log(p*nbin);
        }
        result[v] = fmin(1, fmax(0, kl/log(nbin)));
    }
    mxFree(means);
    mxFree(counts);
}
