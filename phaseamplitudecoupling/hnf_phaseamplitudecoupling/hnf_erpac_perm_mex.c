#include "mex.h"
#include <float.h>
#include <math.h>

/* Across-trial circular-linear correlations. Phase [trial,time] is shared
 * by each high-frequency/channel block of weight [trial,time*blocks].
 * MATLAB supplies one trial order per permutation, shared by all weights.
 * Cache the two-dimensional phase subspace. Keep the original centering
 * order for each permuted linear variable (including rounding behavior). */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3 || nlhs != 1)
        mexErrMsgIdAndTxt("NeuroFreq:ERPACPermArguments", "Use rho = hnf_erpac_perm_mex(phase, weights, orders).");
    for (int k = 0; k < 3; ++k)
        if (!mxIsDouble(prhs[k]) || mxIsComplex(prhs[k]) || mxIsSparse(prhs[k]) ||
            mxIsEmpty(prhs[k]) || mxGetNumberOfDimensions(prhs[k]) != 2)
            mexErrMsgIdAndTxt("NeuroFreq:ERPACPermInput", "Use nonempty full real double matrices.");
    const mwSize n = mxGetM(prhs[0]), times = mxGetN(prhs[0]);
    const mwSize vars = mxGetN(prhs[1]), perms = mxGetN(prhs[2]);
    if (mxGetM(prhs[1]) != n || mxGetM(prhs[2]) != n || vars % times != 0)
        mexErrMsgIdAndTxt("NeuroFreq:ERPACPermSize", "Trial counts must match; weight columns must contain whole time blocks.");
    const double *orders = mxGetPr(prhs[2]), *phase = mxGetPr(prhs[0]), *weights = mxGetPr(prhs[1]);
    mwIndex *order = mxCalloc(mxGetNumberOfElements(prhs[2]),sizeof(*order));
    mwSize *seen = mxCalloc(n,sizeof(*seen));
    for (mwSize p = 0; p < perms; ++p) {
        for (mwSize t = 0; t < n; ++t) {
            double index = orders[t+p*n];
            if (!mxIsFinite(index) || index < 1 || index > (double)n || floor(index) != index)
                mexErrMsgIdAndTxt("NeuroFreq:ERPACPermOrder", "Orders must be one-based trial permutations.");
            mwIndex j = (mwIndex)(index-1);
            if (seen[j] == p+1)
                mexErrMsgIdAndTxt("NeuroFreq:ERPACPermOrder","A trial is repeated in an ERPAC permutation.");
            seen[j] = p+1; order[t+p*n] = j;
        }
    }
    mxFree(seen);
    plhs[0] = mxCreateDoubleMatrix(perms, vars, mxREAL);
    double *out = mxGetPr(plhs[0]);
    for (mwSize i = 0; i < mxGetNumberOfElements(plhs[0]); ++i) out[i] = mxGetNaN();
    if (n < 3) { mxFree(order); return; }
    double *s = mxCalloc(n,sizeof(*s)), *c = mxCalloc(n,sizeof(*c));
    double *w = mxCalloc(n,sizeof(*w)), *x = mxCalloc(n,sizeof(*x));
    for (mwSize time = 0; time < times; ++time) {
        const double *a = phase+time*n;
        int valid = 1;
        double ms = 0, mc = 0, s0 = sin(a[0]), c0 = cos(a[0]);
        for (mwSize t = 0; t < n; ++t) {
            if (!mxIsFinite(a[t])) { valid = 0; break; }
            s[t] = sin(a[t])-s0; c[t] = cos(a[t])-c0; ms += s[t]; mc += c[t];
        }
        if (!valid) continue;
        double ss = 0, cc = 0, sc = 0;
        for (mwSize t = 0; t < n; ++t) {
            s[t] -= ms/n; c[t] -= mc/n;
            ss += s[t]*s[t]; cc += c[t]*c[t]; sc += s[t]*c[t];
        }
        if (fmax(ss,cc) == 0) continue;
        double *first = ss >= cc ? s : c, *other = ss >= cc ? c : s;
        double norm1 = sqrt(fmax(ss,cc)), proj = sc/norm1, residual = 0;
        for (mwSize t = 0; t < n; ++t) {
            first[t] /= norm1; other[t] -= proj*first[t]; residual += other[t]*other[t];
        }
        double norm2 = sqrt(residual);
        if (norm2 <= n*DBL_EPSILON*norm1) continue;

        for (mwSize v = time; v < vars; v += times) {
            const double *y = weights+v*n;
            double scale = 0;
            for (mwSize t = 0; t < n; ++t) {
                if (!mxIsFinite(y[t])) { scale = 0; break; }
                if (fabs(y[t]) > scale) scale = fabs(y[t]);
            }
            if (scale == 0) continue;
            for (mwSize t = 0; t < n; ++t) w[t] = y[t]/scale;
            for (mwSize p = 0; p < perms; ++p) {
                const mwIndex *indices = order+p*n;
                double mean = 0, xx = 0, du = 0, dv = 0;
                for (mwSize t = 0; t < n; ++t) { x[t] = w[indices[t]]-w[indices[0]]; mean += x[t]; }
                for (mwSize t = 0; t < n; ++t) {
                    x[t] -= mean/n; xx += x[t]*x[t];
                    du += x[t]*first[t]; dv += x[t]*other[t];
                }
                if (xx > 0) out[p+v*perms] = fmin(1,hypot(du,dv/norm2)/sqrt(xx));
            }
        }
    }
    mxFree(order); mxFree(s); mxFree(c); mxFree(x); mxFree(w);
}
