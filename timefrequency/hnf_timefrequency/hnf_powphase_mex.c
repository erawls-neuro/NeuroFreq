/* Squared magnitude and phase in one pass, used after TF convolutions.
 * Power is computed in double even for single Fourier coefficients. */
#include "mex.h"
#include <math.h>
#include <float.h>

static double power_value(double x,double y)
{
    double m = fmax(fabs(x),fabs(y));
    if ((m!=0.0 && m<sqrt(DBL_MIN)) || !isfinite(x) || !isfinite(y)) {
        double mag = hypot(x,y);
        return mag*mag;
    }
    return x*x+y*y;
}

void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *in[])
{
    if (nlhs!=2 || nrhs!=1 || (!mxIsSingle(in[0]) && !mxIsDouble(in[0])) || mxIsSparse(in[0]))
        mexErrMsgIdAndTxt("NeuroFreq:PowPhase","Expected one full floating-point array and two outputs.");
    const mxArray *a=in[0];
    out[0]=mxCreateNumericArray(mxGetNumberOfDimensions(a),mxGetDimensions(a),mxDOUBLE_CLASS,mxREAL);
    out[1]=mxCreateNumericArray(mxGetNumberOfDimensions(a),mxGetDimensions(a),mxGetClassID(a),mxREAL);
    double *p=mxGetPr(out[0]); size_t n=mxGetNumberOfElements(a);
    if (mxIsSingle(a)) {
        float *phase=(float *)mxGetData(out[1]);
#if MX_HAS_INTERLEAVED_COMPLEX
        const mxComplexSingle *z=mxIsComplex(a)?mxGetComplexSingles(a):NULL;
        const float *r=z?NULL:mxGetSingles(a);
#else
        const float *r=(const float *)mxGetData(a),*im=(const float *)mxGetImagData(a);
#endif
        for (size_t k=0;k<n;k++) {
#if MX_HAS_INTERLEAVED_COMPLEX
            float x=z?z[k].real:r[k],y=z?z[k].imag:0;
#else
            float x=r[k],y=im?im[k]:0;
#endif
            /* Every finite single-precision square fits in double. */
            p[k]=isfinite(x) && isfinite(y) ? (double)x*x+(double)y*y : power_value(x,y);
            phase[k]=atan2f(y,x);
        }
    } else {
        double *phase=mxGetPr(out[1]);
#if MX_HAS_INTERLEAVED_COMPLEX
        const mxComplexDouble *z=mxIsComplex(a)?mxGetComplexDoubles(a):NULL;
        const double *r=z?NULL:mxGetDoubles(a);
#else
        const double *r=mxGetPr(a),*im=mxGetPi(a);
#endif
        for (size_t k=0;k<n;k++) {
#if MX_HAS_INTERLEAVED_COMPLEX
            double x=z?z[k].real:r[k],y=z?z[k].imag:0;
#else
            double x=r[k],y=im?im[k]:0;
#endif
            p[k]=power_value(x,y);phase[k]=atan2(y,x);
        }
    }
}
