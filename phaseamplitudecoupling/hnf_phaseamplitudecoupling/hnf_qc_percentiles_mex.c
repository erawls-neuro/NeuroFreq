/* NF's channel-wise PSD envelopes. Each frequency column is sorted once for
 * the three reported quantiles. Exclude nonfinite channels independently
 * at each frequency; NF uses linear interpolation at 1 + (n-1)*p/100. */
#include "mex.h"
#include <stdlib.h>

static int ascending(const void *left, const void *right)
{
    const double a = *(const double *)left, b = *(const double *)right;
    return (a > b) - (a < b);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const double percentages[] = {50, 10, 90};
    const double *data;
    double *out, *work;
    mwSize channels, frequencies, freq, channel, k;
    if (nlhs != 1 || nrhs != 1 || !mxIsDouble(prhs[0]) || mxIsSparse(prhs[0]) ||
        mxIsComplex(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgIdAndTxt("NeuroFreq:QCPercentiles", "Expected a real full channels-by-frequencies matrix.");
    channels = mxGetM(prhs[0]);
    frequencies = mxGetN(prhs[0]);
    data = mxGetPr(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(3, frequencies, mxREAL);
    out = mxGetPr(plhs[0]);
    work = mxCalloc(channels, sizeof(double));
    for (freq = 0; freq < frequencies; ++freq) {
        mwSize count = 0;
        for (channel = 0; channel < channels; ++channel) {
            double value = data[freq * channels + channel];
            if (mxIsFinite(value)) work[count++] = value;
        }
        if (count > 1) qsort(work, count, sizeof(double), ascending);
        for (k = 0; k < 3; ++k) {
            if (count < 2) {
                out[3 * freq + k] = count ? work[0] : mxGetNaN();
            } else {
                double position = 1 + (count - 1) * percentages[k] / 100;
                mwIndex low = (mwIndex)position - 1;
                double fraction = position - (low + 1);
                mwIndex high = low + (fraction > 0);
                const double a = work[low], b = work[high];
                const double delta = b - a;
                if (!mxIsFinite(delta))
                    out[3 * freq + k] = (1.0 - fraction) * a + fraction * b;
                else if (fraction > 0.5)
                    out[3 * freq + k] = b - (1.0 - fraction) * delta;
                else
                    out[3 * freq + k] = a + fraction * delta;
            }
        }
    }
    mxFree(work);
}
