/* Finite-scale norm/RMS along any array dimension. The wrapper retains
 * complex inputs and diagnostics; real CPU arrays can avoid array-sized
 * masks, normalized copies, and squared-amplitude temporaries here. */
#include "mex.h"
#include <math.h>

template<typename T>
static void reduce(const T *in, T *out, mwSize outSize, mwSize stride,
                   mwSize axisSize, bool omit, bool rmsMode)
{
    for (mwSize k = 0; k < outSize; ++k) {
        mwSize base = (k / stride) * axisSize * stride + k % stride;
        mwSize count = 0;
        double scale = 0.0;
        bool hasInf = false, hasNan = false;
        for (mwSize j = 0; j < axisSize; ++j) {
            double a = fabs((double)in[base + j*stride]);
            if (isnan(a)) {
                hasNan = true;
                continue;
            }
            ++count;
            if (isinf(a)) hasInf = true;
            if (a > scale) scale = a;
        }
        if ((hasNan && !omit) || (rmsMode && count == 0)) {
            out[k] = (T)mxGetNaN();
        } else if (hasInf) {
            out[k] = (T)mxGetInf();
        } else if (scale == 0.0) {
            out[k] = 0;
        } else {
            double squares = 0.0;
            for (mwSize j = 0; j < axisSize; ++j) {
                double a = (double)in[base + j*stride];
                if (!isnan(a)) {
                    double scaled = a / scale;
                    squares += scaled * scaled;
                }
            }
            out[k] = (T)(scale * sqrt(rmsMode ? squares / (double)count : squares));
        }
    }
}

void mexFunction(int nlhs,mxArray **plhs,int nrhs,const mxArray **prhs)
{
    const mxArray *x;
    const mwSize *dims;
    mwSize *shape, nd, k, axis, axisSize=1, stride=1, outSize;
    double dim;
    int omit, rmsMode;
    if (nrhs!=4 || nlhs!=1)
        mexErrMsgIdAndTxt("NeuroFreq:length:arity","Expected X, dimension and two logical flags.");
    x=prhs[0];
    if (!(mxIsDouble(x)||mxIsSingle(x)) || mxIsSparse(x)||mxIsComplex(x) ||
        !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetNumberOfElements(prhs[1])!=1 ||
        !mxIsLogicalScalar(prhs[2]) || !mxIsLogicalScalar(prhs[3]))
        mexErrMsgIdAndTxt("NeuroFreq:length:type","Expected a full real floating-point array, dimension and scalar flags.");
    dim=mxGetScalar(prhs[1]);
    if (!isfinite(dim) || dim<1.0 || floor(dim)!=dim)
        mexErrMsgIdAndTxt("NeuroFreq:length:dim","Dimension must be a positive integer.");
    omit=mxIsLogicalScalarTrue(prhs[2]); rmsMode=mxIsLogicalScalarTrue(prhs[3]);
    dims=mxGetDimensions(x); nd=mxGetNumberOfDimensions(x);
    /* An axis past nd is a singleton and preserves the input shape. */
    axis = dim>(double)nd ? nd : (mwSize)dim-1;
    shape=(mwSize *)mxMalloc(nd*sizeof(mwSize));
    for (k=0; k<nd; ++k) {
        shape[k]=dims[k];
        if (k<axis) stride*=dims[k];
        if (k==axis) {axisSize=dims[k]; shape[k]=1;}
    }
    plhs[0]=mxCreateNumericArray(nd,shape,mxGetClassID(x),mxREAL); mxFree(shape);
    outSize=mxGetNumberOfElements(plhs[0]);
    if (mxIsDouble(x)) {
        reduce(mxGetPr(x), mxGetPr(plhs[0]), outSize, stride, axisSize, omit, rmsMode);
    } else {
        reduce((const float *)mxGetData(x), (float *)mxGetData(plhs[0]),
               outSize, stride, axisSize, omit, rmsMode);
    }
}
