#include "mex.h"
#include <math.h>
#include <stdint.h>
#include <string.h>

/* Native equivalent of hnf_zero_phase.m. Independent columns are traces.
 * Reuse one reflected trace rather than allocating a matrix per section.
 * Reverse the traversal direction in-place instead of flipping a matrix.
 * The caller has already distributed the Butterworth gain over the cascade. */
typedef struct {
    double b0, b1, b2, a1, a2, initial1, initial2;
} Section;

static void filter_trace(double *trace, const size_t length,
                         const Section *sections, const size_t count)
{
    for (int reverse = 0; reverse < 2; ++reverse) {
        for (size_t s = 0; s < count; ++s) {
            Section c = sections[s];
            double start = trace[reverse ? length-1 : 0];
            double z1 = c.initial1 * start, z2 = c.initial2 * start;
            for (size_t offset = 0; offset < length; ++offset) {
                size_t t = reverse ? length-1-offset : offset;
                double x = trace[t];
                double y = c.b0*x + z1;
                z1 = (z2 + c.b1*x) - c.a1*y;
                z2 = c.b2*x - c.a2*y;
                trace[t] = y;
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3 || nlhs != 1)
        mexErrMsgIdAndTxt("NeuroFreq:SOSArgs","Use y = hnf_sosfiltfilt_mex(sos,pad,x).");
    const mxArray *sos = prhs[0], *p = prhs[1], *x = prhs[2];
    if (!mxIsDouble(sos) || mxIsSparse(sos) || mxIsComplex(sos) ||
        mxGetNumberOfDimensions(sos) != 2 || !mxGetM(sos) || mxGetN(sos) != 6 ||
        !mxIsDouble(x) || mxIsSparse(x) || mxGetNumberOfDimensions(x) != 2 ||
        mxIsEmpty(x) || !mxIsDouble(p) || mxIsSparse(p) || mxIsComplex(p) ||
        mxGetNumberOfElements(p) != 1)
        mexErrMsgIdAndTxt("NeuroFreq:SOSInput","Expected a real S-by-6 coefficient matrix, integer pad and full double input.");
    size_t n = mxGetM(x), columns = mxGetN(x), nSections = mxGetM(sos);
    double pvalue = mxGetScalar(p);
    if (!(pvalue >= 0 && pvalue < (double)n && pvalue == floor(pvalue)))
        mexErrMsgIdAndTxt("NeuroFreq:SOSPad","The integer padding must be shorter than an epoch.");
    size_t pad = (size_t)pvalue;
    if (n > SIZE_MAX/sizeof(double) || pad > (SIZE_MAX/sizeof(double)-n)/2)
        mexErrMsgIdAndTxt("NeuroFreq:SOSSize","The reflected trace is too large.");
    size_t length = n + 2*pad;
    const double *coefficients = mxGetPr(sos);
    Section *sections = mxCalloc(nSections,sizeof(*sections));
    for (size_t s = 0; s < nSections; ++s) {
        double a0 = coefficients[3*nSections+s];
        if (a0 == 0)
            mexErrMsgIdAndTxt("NeuroFreq:SOSCoefficient","Invalid denominator coefficient.");
        Section c;
        c.b0 = coefficients[s]/a0;
        c.b1 = coefficients[nSections+s]/a0;
        c.b2 = coefficients[2*nSections+s]/a0;
        c.a1 = coefficients[4*nSections+s]/a0;
        c.a2 = coefficients[5*nSections+s]/a0;
        double dc = (c.b0+c.b1+c.b2)/(1+c.a1+c.a2);
        c.initial1 = dc-c.b0;
        c.initial2 = c.b2-c.a2*dc;
        sections[s] = c;
    }
    double *trace = mxMalloc(length*sizeof(*trace));
    int complexX = mxIsComplex(x);
    plhs[0] = mxCreateDoubleMatrix(n,columns,complexX ? mxCOMPLEX : mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
    const mxComplexDouble *z = complexX ? mxGetComplexDoubles(x) : NULL;
    mxComplexDouble *outZ = complexX ? mxGetComplexDoubles(plhs[0]) : NULL;
    const double *xr = complexX ? NULL : mxGetDoubles(x);
    double *yr = complexX ? NULL : mxGetDoubles(plhs[0]);
#else
    const double *xr = mxGetPr(x), *xi = mxGetPi(x);
    double *yr = mxGetPr(plhs[0]), *yi = mxGetPi(plhs[0]);
#endif
    for (size_t col = 0; col < columns; ++col) {
        size_t base = col*n;
        for (int imag = 0; imag <= complexX; ++imag) {
            double *center = trace+pad;
#if MX_HAS_INTERLEAVED_COMPLEX
            if (complexX) {
                for (size_t t = 0; t < n; ++t)
                    center[t] = imag ? z[base+t].imag : z[base+t].real;
            } else memcpy(center,xr+base,n*sizeof(*center));
#else
            memcpy(center,(imag ? xi : xr)+base,n*sizeof(*center));
#endif
            for (size_t t = 0; t < pad; ++t) {
                trace[t] = 2*center[0]-center[pad-t];
                center[n+t] = 2*center[n-1]-center[n-2-t];
            }
            filter_trace(trace,length,sections,nSections);
#if MX_HAS_INTERLEAVED_COMPLEX
            if (complexX) {
                for (size_t t = 0; t < n; ++t) {
                    if (imag) outZ[base+t].imag = center[t];
                    else outZ[base+t].real = center[t];
                }
            } else memcpy(yr+base,center,n*sizeof(*center));
#else
            memcpy((imag ? yi : yr)+base,center,n*sizeof(*center));
#endif
        }
    }
    mxFree(trace);
    mxFree(sections);
}
