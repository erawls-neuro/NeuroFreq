#include "mex.h"
#include <math.h>

/*
 * hnf_mvl_mex
 *
 * mvl = hnf_mvl_mex(phas, pow)
 *
 * Mean Vector Length (MVL) kernel:
 *   - First dimension = samples (time)
 *   - All trailing dimensions define independent phase–amplitude series
 *   - Output has same size as inputs, but with first dim collapsed to 1
 *
 * MVL = | mean_s( pow(s) * exp(1i * phas(s)) ) |
 *      = sqrt( (sum_s pow(s)*cos(phas(s)))^2 + (sum_s pow(s)*sin(phas(s)))^2 ) / nSamples
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Inputs */
    const mxArray *phas_mx;
    const mxArray *pow_mx;

    double *phas;
    double *pow;

    /* Size info */
    mwSize nDims;
    const mwSize *dims;
    mwSize i;
    mwSize nSamples;
    mwSize nVars;
    mwSize totalElems;

    /* Output */
    mxArray *mvl_mx;
    double *mvl_out;
    mwSize *dimsOut;

    /* Indices and accumulators */
    mwSize v;
    mwSize s;

    /* --- Argument checks --- */

    if (nrhs != 2)
    {
        mexErrMsgIdAndTxt("hnf_mvl_mex:nrhs",
                          "Two input arguments required: phas, pow.");
    }

    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("hnf_mvl_mex:nlhs",
                          "Too many output arguments.");
    }

    phas_mx = prhs[0];
    pow_mx  = prhs[1];

    if (!mxIsDouble(phas_mx) || mxIsComplex(phas_mx))
    {
        mexErrMsgIdAndTxt("hnf_mvl_mex:phasType",
                          "phas must be real double.");
    }

    if (!mxIsDouble(pow_mx) || mxIsComplex(pow_mx))
    {
        mexErrMsgIdAndTxt("hnf_mvl_mex:powType",
                          "pow must be real double.");
    }

    nDims = mxGetNumberOfDimensions(phas_mx);
    dims  = mxGetDimensions(phas_mx);

    if (nDims < 1)
    {
        mexErrMsgIdAndTxt("hnf_mvl_mex:dimensions",
                          "phas must have at least one dimension.");
    }

    {
        mwSize nDimsPow = mxGetNumberOfDimensions(pow_mx);
        const mwSize *dimsPow = mxGetDimensions(pow_mx);

        if (nDimsPow != nDims)
        {
            mexErrMsgIdAndTxt("hnf_mvl_mex:sizeMismatch",
                              "phas and pow must have the same size.");
        }

        for (i = 0; i < nDims; ++i)
        {
            if (dims[i] != dimsPow[i])
            {
                mexErrMsgIdAndTxt("hnf_mvl_mex:sizeMismatch",
                                  "phas and pow must have the same size.");
            }
        }
    }

    /* Size info derived from dims */
    nSamples   = dims[0];
    totalElems = mxGetNumberOfElements(phas_mx);

    if (nSamples < 1)
    {
        mexErrMsgIdAndTxt("hnf_mvl_mex:nSamples",
                          "first dimension must have at least one sample.");
    }

    if (totalElems % nSamples != 0)
    {
        mexErrMsgIdAndTxt("hnf_mvl_mex:internal",
                          "Number of elements not divisible by nSamples.");
    }

    nVars = totalElems / nSamples;

    /* Pointers to input data */
    phas = mxGetPr(phas_mx);
    pow  = mxGetPr(pow_mx);

    /* --- Create output array with first dim collapsed to 1 --- */

    dimsOut = (mwSize *)mxCalloc((size_t)nDims, sizeof(mwSize));
    if (dimsOut == NULL)
    {
        mexErrMsgIdAndTxt("hnf_mvl_mex:alloc",
                          "Output dimension allocation failed.");
    }

    dimsOut[0] = 1;
    for (i = 1; i < nDims; ++i)
    {
        dimsOut[i] = dims[i];
    }

    mvl_mx = mxCreateNumericArray(nDims, dimsOut, mxDOUBLE_CLASS, mxREAL);
    if (mvl_mx == NULL)
    {
        mexErrMsgIdAndTxt("hnf_mvl_mex:alloc",
                          "Output array allocation failed.");
    }

    mvl_out = mxGetPr(mvl_mx);

    /* --- Main MVL computation per variable --- */

    for (v = 0; v < nVars; ++v)
    {
        mwSize base = v * nSamples;

        double sumRe = 0.0;
        double sumIm = 0.0;

        for (s = 0; s < nSamples; ++s)
        {
            mwSize idx = base + s;

            double phi    = phas[idx];
            double ampl   = pow[idx];
            double cphi   = cos(phi);
            double sphi   = sin(phi);

            sumRe += ampl * cphi;
            sumIm += ampl * sphi;
        }

        if (nSamples > 0)
        {
            double meanRe = sumRe / (double)nSamples;
            double meanIm = sumIm / (double)nSamples;
            double mag    = sqrt(meanRe * meanRe + meanIm * meanIm);

            mvl_out[v] = mag;
        }
        else
        {
            mvl_out[v] = mxGetNaN();
        }
    }

    plhs[0] = mvl_mx;

    mxFree(dimsOut);
}