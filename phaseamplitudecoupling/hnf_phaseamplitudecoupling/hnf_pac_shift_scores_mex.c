#include "mex.h"
#include <math.h>
#include <string.h>

/* Score the specified positive circular shifts. MATLAB owns random sampling.
 * Columns are independent channel/trial series. The metric's phase geometry
 * and the weight scale are invariant to shifting. Compute them once and
 * visit samples in the same order as hnf_mi_mex / hnf_mvl_mex. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4 || nlhs != 1)
        mexErrMsgIdAndTxt("NeuroFreq:ShiftScoresArguments",
            "Use scores = hnf_pac_shift_scores_mex(phase, weights, shifts, nbin).");
    for (int k = 0; k < nrhs; ++k)
        if (!mxIsDouble(prhs[k]) || mxIsComplex(prhs[k]) || mxIsSparse(prhs[k]) ||
            mxIsEmpty(prhs[k]) || mxGetNumberOfDimensions(prhs[k]) != 2)
            mexErrMsgIdAndTxt("NeuroFreq:ShiftScoresInput", "Use nonempty full real double matrices.");
    const mwSize n = mxGetM(prhs[0]), vars = mxGetN(prhs[0]), perms = mxGetN(prhs[2]);
    if (mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != vars || mxGetM(prhs[2]) != vars ||
        mxGetNumberOfElements(prhs[3]) != 1)
        mexErrMsgIdAndTxt("NeuroFreq:ShiftScoresSize", "Weights must match phase; shifts must be series by permutations.");
    const double nbin = mxGetScalar(prhs[3]);
    if (!mxIsFinite(nbin) || nbin < 0 || nbin == 1 || floor(nbin) != nbin)
        mexErrMsgIdAndTxt("NeuroFreq:ShiftScoresBins", "Use zero for MVL or an integer of at least two for MI.");
    const double *shift = mxGetPr(prhs[2]);
    for (mwSize k = 0; k < mxGetNumberOfElements(prhs[2]); ++k)
        if (!mxIsFinite(shift[k]) || shift[k] < 0 || shift[k] >= (double)n || floor(shift[k]) != shift[k])
            mexErrMsgIdAndTxt("NeuroFreq:ShiftScoresShift", "Shifts must be integers from zero to the sample count minus one.");

    plhs[0] = mxCreateDoubleMatrix(perms, vars, mxREAL);
    double *out = mxGetPr(plhs[0]);
    for (mwSize k = 0; k < mxGetNumberOfElements(plhs[0]); ++k) out[k] = mxGetNaN();
    if (nbin > (double)n) return; /* MI with a necessarily empty bin. */
    const mwSize bins = (mwSize)nbin;
    const double pi = 3.14159265358979323846, period = 2*pi;
    const double logBins = bins ? log(nbin) : 0;
    double *w = mxCalloc(n, sizeof(*w));
    mwIndex *bin = bins ? mxCalloc(n, sizeof(*bin)) : NULL;
    mwSize *counts = bins ? mxCalloc(bins, sizeof(*counts)) : NULL;
    double *means = bins ? mxCalloc(bins, sizeof(*means)) : NULL;
    double *cosine = bins ? NULL : mxCalloc(n, sizeof(*cosine));
    double *sine = bins ? NULL : mxCalloc(n, sizeof(*sine));
    const double *phase = mxGetPr(prhs[0]), *weight = mxGetPr(prhs[1]);

    for (mwSize v = 0; v < vars; ++v) {
        const double *a = phase + v*n, *x = weight + v*n;
        double scale = 0;
        int valid = 1;
        for (mwSize t = 0; t < n; ++t) {
            if (!mxIsFinite(a[t]) || !mxIsFinite(x[t]) || (bins && x[t] < 0)) { valid = 0; break; }
            if (fabs(x[t]) > scale) scale = fabs(x[t]);
        }
        if (!valid || (bins && scale == 0)) continue;
        if (scale == 0) scale = 1;
        if (bins) memset(counts, 0, bins*sizeof(*counts));
        for (mwSize t = 0; t < n; ++t) {
            w[t] = x[t]/scale;
            if (bins) {
                double angle = fmod(a[t]+pi, period);
                if (angle < 0) angle += period;
                mwIndex b = (mwIndex)floor(angle/(period/nbin));
                bin[t] = b < bins ? b : bins-1;
                ++counts[bin[t]];
            } else {
                cosine[t] = cos(a[t]); sine[t] = sin(a[t]);
            }
        }
        for (mwSize b = 0; b < bins; ++b) if (counts[b] == 0) valid = 0;
        if (!valid) continue;
        for (mwSize p = 0; p < perms; ++p) {
            const mwSize lag = (mwSize)shift[v+p*vars];
            mwSize j = lag ? n-lag : 0; /* w[(t-lag) mod n]. */
            if (bins) {
                memset(means, 0, bins*sizeof(*means));
                for (mwSize t = 0; t < n; ++t) {
                    means[bin[t]] += w[j];
                    if (++j == n) j = 0;
                }
                double total = 0, kl = 0;
                for (mwSize b = 0; b < bins; ++b) { means[b] /= (double)counts[b]; total += means[b]; }
                if (total == 0) continue;
                for (mwSize b = 0; b < bins; ++b) {
                    double probability = means[b]/total;
                    if (probability > 0) kl += probability*log(probability*nbin);
                }
                out[p+v*perms] = fmin(1, fmax(0, kl/logBins));
            } else {
                double re = 0, im = 0;
                for (mwSize t = 0; t < n; ++t) {
                    re += w[j]*cosine[t]; im += w[j]*sine[t];
                    if (++j == n) j = 0;
                }
                out[p+v*perms] = scale*hypot(re/n, im/n);
            }
        }
    }
    mxFree(w);
    if (bins) { mxFree(bin); mxFree(counts); mxFree(means); }
    else { mxFree(cosine); mxFree(sine); }
}
