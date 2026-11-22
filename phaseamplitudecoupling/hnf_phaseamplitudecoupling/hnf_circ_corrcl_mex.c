#include "mex.h"
#include <float.h>
#include <math.h>

/* See hnf_circ_corrcl_vec.m for the orthogonal projection estimator. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2 || nlhs > 2)
        mexErrMsgIdAndTxt("NeuroFreq:CircArgs","Use [rho,p] = hnf_circ_corrcl_mex(phase,x).");
    for (int k = 0; k < 2; ++k)
        if (!mxIsDouble(prhs[k]) || mxIsComplex(prhs[k]) || mxIsSparse(prhs[k]))
            mexErrMsgIdAndTxt("NeuroFreq:CircInput","Both inputs must be full, real double arrays.");
    mwSize nd = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *dims = mxGetDimensions(prhs[0]);
    if (mxIsEmpty(prhs[0]) || mxGetNumberOfDimensions(prhs[1]) != nd)
        mexErrMsgIdAndTxt("NeuroFreq:CircInput","Nonempty arrays of equal size required.");
    for (mwSize k = 0; k < nd; ++k)
        if (dims[k] != mxGetDimensions(prhs[1])[k])
            mexErrMsgIdAndTxt("NeuroFreq:CircInput","Phase and x must have equal size.");
    mwSize n = dims[0], nvars = mxGetNumberOfElements(prhs[0])/n;
    mwSize *shape = mxCalloc(nd, sizeof(*shape));
    for (mwSize k = 0; k < nd; ++k) shape[k] = k ? dims[k] : 1;
    plhs[0] = mxCreateNumericArray(nd,shape,mxDOUBLE_CLASS,mxREAL);
    if (nlhs > 1) plhs[1] = mxCreateNumericArray(nd,shape,mxDOUBLE_CLASS,mxREAL);
    mxFree(shape);
    double *rho = mxGetPr(plhs[0]), *pval = nlhs > 1 ? mxGetPr(plhs[1]) : NULL;
    const double *phase = mxGetPr(prhs[0]), *linear = mxGetPr(prhs[1]);
    double *s = mxCalloc(n,sizeof(*s)), *c = mxCalloc(n,sizeof(*c)), *x = mxCalloc(n,sizeof(*x));
    for (mwSize v = 0; v < nvars; ++v) {
        rho[v] = mxGetNaN();
        if (pval) pval[v] = mxGetNaN();
        if (n < 3) continue;
        const double *a = phase+v*n, *y = linear+v*n;
        double scale = 0, ms = 0, mc = 0, mx = 0;
        for (mwSize t = 0; t < n; ++t) {
            if (!mxIsFinite(a[t]) || !mxIsFinite(y[t])) { scale = 0; break; }
            if (fabs(y[t]) > scale) scale = fabs(y[t]);
        }
        if (scale == 0) continue;
        double s0 = sin(a[0]), c0 = cos(a[0]), x0 = y[0]/scale;
        for (mwSize t = 0; t < n; ++t) {
            s[t] = sin(a[t])-s0; c[t] = cos(a[t])-c0; x[t] = y[t]/scale-x0;
            ms += s[t]; mc += c[t]; mx += x[t];
        }
        double ss = 0, cc = 0, xx = 0, sc = 0;
        for (mwSize t = 0; t < n; ++t) {
            s[t] -= ms/n; c[t] -= mc/n; x[t] -= mx/n;
            ss += s[t]*s[t]; cc += c[t]*c[t]; xx += x[t]*x[t]; sc += s[t]*c[t];
        }
        if (xx == 0 || fmax(ss,cc) == 0) continue;
        double *first = ss >= cc ? s : c, *other = ss >= cc ? c : s;
        double norm1 = sqrt(fmax(ss,cc)), proj = sc/norm1, residualNorm2 = 0;
        double du = 0, dv = 0;
        for (mwSize t = 0; t < n; ++t) {
            first[t] /= norm1; other[t] -= proj*first[t];
            residualNorm2 += other[t]*other[t];
            du += x[t]*first[t]; dv += x[t]*other[t];
        }
        double norm2 = sqrt(residualNorm2);
        if (norm2 <= n*DBL_EPSILON*norm1) continue;
        rho[v] = fmin(1,hypot(du,dv/norm2)/sqrt(xx));
        if (pval) pval[v] = exp(-0.5*n*rho[v]*rho[v]);
    }
    mxFree(s); mxFree(c); mxFree(x);
}
