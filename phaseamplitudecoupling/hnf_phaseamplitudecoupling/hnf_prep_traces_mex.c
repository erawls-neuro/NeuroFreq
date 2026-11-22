/* Prepare one epoch (channels x samples). MATLAB retains its Hilbert FFT.
 * Every channel shares the same sample grid. On that grid, 1, u and
 * u^2 - mean(u^2) are mutually orthogonal for u equally spaced in [-1,1].
 * Projecting onto these columns removes the same quadratic as polyfit
 * without factoring a new Vandermonde matrix for every channel/trial.
 * Accumulate centered values to retain precision in signals with a DC offset. */
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const double *taper;
    double *u, *q, *result;
    const double *input;
    mwSize channels, samples, channel, time;
    double meanSquare = 0, normLinear = 0, normQuad = 0;
    if (nrhs != 2 || nlhs != 1 || !mxIsDouble(prhs[0]) ||
        mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) ||
        mxGetNumberOfDimensions(prhs[0]) != 2 ||
        !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]))
        mexErrMsgIdAndTxt("NeuroFreq:PrepArguments", "Expected a centered real double epoch and a double taper or [].");
    channels = mxGetM(prhs[0]);
    samples = mxGetN(prhs[0]);
    if (samples < 3 || (mxGetNumberOfElements(prhs[1]) != 0 &&
                        (mwSize)mxGetNumberOfElements(prhs[1]) != samples))
        mexErrMsgIdAndTxt("NeuroFreq:PrepGrid", "Expected at least three samples and one taper weight per sample.");
    input = mxGetPr(prhs[0]);
    taper = mxGetNumberOfElements(prhs[1]) ? mxGetPr(prhs[1]) : NULL;
    plhs[0] = mxCreateDoubleMatrix(channels, samples, mxREAL);
    result = mxGetPr(plhs[0]);
    u = mxCalloc(samples, sizeof(double));
    q = mxCalloc(samples, sizeof(double));
    for (time = 0; time < samples; ++time) {
        u[time] = (2.0 * time - (samples - 1)) / (samples - 1);
        normLinear += u[time] * u[time];
    }
    meanSquare = normLinear / samples;
    for (time = 0; time < samples; ++time) {
        q[time] = u[time] * u[time] - meanSquare;
        normQuad += q[time] * q[time];
    }
    for (channel = 0; channel < channels; ++channel) {
        double offset = 0, linear = 0, quadratic = 0;
        for (time = 0; time < samples; ++time)
            offset += input[channel + time * channels];
        offset /= samples;
        for (time = 0; time < samples; ++time) {
            double y = input[channel + time * channels] - offset;
            linear += y * u[time];
            quadratic += y * q[time];
        }
        linear /= normLinear;
        quadratic /= normQuad;
        for (time = 0; time < samples; ++time) {
            mwSize index = channel + time * channels;
            double y = input[index] - offset;
            y -= linear * u[time] + quadratic * q[time];
            result[index] = taper ? y * taper[time] : y;
        }
    }
    mxFree(u); mxFree(q);
}
