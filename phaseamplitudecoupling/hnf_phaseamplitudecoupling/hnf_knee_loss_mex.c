/* Log-domain knee objective for repeated Nelder-Mead evaluations.
 * Parameters: offset, knee and exponent. Spectra are in log10 power.
 * The exact same loss is available in hnf_error_expo_function.m. */
#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nlhs != 1 || nrhs != 3)
        mexErrMsgIdAndTxt("NeuroFreq:knee:arity","Expected parameters, log frequencies and spectrum.");
    for (int i = 0; i < 3; ++i)
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]))
            mexErrMsgIdAndTxt("NeuroFreq:knee:type","Expected dense real doubles.");
    size_t bins = mxGetNumberOfElements(prhs[1]);
    if (mxGetNumberOfElements(prhs[0]) != 3 || mxGetNumberOfElements(prhs[2]) != bins)
        mexErrMsgIdAndTxt("NeuroFreq:knee:size","Expected three parameters and matching observation counts.");
    const double *p = mxGetPr(prhs[0]), *logFreq = mxGetPr(prhs[1]), *y = mxGetPr(prhs[2]);
    double logKnee = log(fabs(p[1])), loss = 0;
    const double logTen = log(10.0);
    for (size_t i = 0; i < bins; ++i) {
        double logPower = p[2]*logFreq[i], logSum = logPower;
        if (p[1] != 0)
            logSum = fmax(logKnee,logPower) + log1p(exp(-fabs(logKnee-logPower)));
        double r = y[i] - (p[0]-logSum/logTen);
        loss += r*r;
    }
    plhs[0] = mxCreateDoubleScalar(loss);
}
