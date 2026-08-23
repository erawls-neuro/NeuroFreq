#include "mex.h"
#include <math.h>

/*
 * hnf_mi_mex
 *
 * mi = hnf_mi_mex(phas, pow, nbin)
 *
 * Vectorized Tort MI kernel:
 *  - First dim = samples (time)
 *  - All trailing dims define independent phase–amplitude pairs
 *  - Output has the same size as inputs, but with first dim collapsed to 1
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Inputs */
    const mxArray *phas_mx;
    const mxArray *pow_mx;
    const mxArray *nbin_mx;

    double *phas;
    double *pow;
    double nbin_d;
    mwSize nbin_i;

    /* Size info */
    mwSize nDims;
    const mwSize *dims;
    mwSize i;
    mwSize nSamples;
    mwSize nVars;
    mwSize totalElems;

    /* Accumulators */
    double *sumAmp;
    mwSize *countAmp;
    double *colSum;
    double *distKL;
    double *mi_vec;

    /* Output */
    mxArray *mi_mx;
    double *mi_out;
    mwSize *dimsOut;

    /* Constants */
    const double PI = 3.14159265358979323846;
    double binWidth;
    double log_nbin;
    double NaNval;

    /* Indices */
    mwSize v;
    mwSize s;
    mwSize b;

    /* --- Argument checks --- */

    if (nrhs != 3)
    {
        mexErrMsgIdAndTxt("hnf_mi_mex:nrhs",
                          "Three input arguments required: phas, pow, nbin.");
    }
    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("hnf_mi_mex:nlhs",
                          "Too many output arguments.");
    }

    phas_mx = prhs[0];
    pow_mx  = prhs[1];
    nbin_mx = prhs[2];

    if (!mxIsDouble(phas_mx) || mxIsComplex(phas_mx))
    {
        mexErrMsgIdAndTxt("hnf_mi_mex:phasType",
                          "phas must be real double.");
    }
    if (!mxIsDouble(pow_mx) || mxIsComplex(pow_mx))
    {
        mexErrMsgIdAndTxt("hnf_mi_mex:powType",
                          "pow must be real double.");
    }

    nDims = mxGetNumberOfDimensions(phas_mx);
    dims  = mxGetDimensions(phas_mx);

    if (nDims < 1)
    {
        mexErrMsgIdAndTxt("hnf_mi_mex:dimensions",
                          "phas must have at least one dimension.");
    }

    {
        mwSize nDimsPow = mxGetNumberOfDimensions(pow_mx);
        const mwSize *dimsPow = mxGetDimensions(pow_mx);
        if (nDimsPow != nDims)
        {
            mexErrMsgIdAndTxt("hnf_mi_mex:sizeMismatch",
                              "phas and pow must have the same size.");
        }
        for (i = 0; i < nDims; ++i)
        {
            if (dims[i] != dimsPow[i])
            {
                mexErrMsgIdAndTxt("hnf_mi_mex:sizeMismatch",
                                  "phas and pow must have the same size.");
            }
        }
    }

    /* nbin checks */
    if (!mxIsDouble(nbin_mx) || mxIsComplex(nbin_mx) || mxGetNumberOfElements(nbin_mx) != 1)
    {
        mexErrMsgIdAndTxt("hnf_mi_mex:nbinType",
                          "nbin must be a real scalar.");
    }

    nbin_d = mxGetScalar(nbin_mx);

    if (nbin_d < 1.0 || floor(nbin_d) != nbin_d)
    {
        mexErrMsgIdAndTxt("hnf_mi_mex:nbinValue",
                          "nbin must be a positive integer scalar.");
    }

    nbin_i = (mwSize)nbin_d;

    /* Size info */
    nSamples   = dims[0];
    totalElems = mxGetNumberOfElements(phas_mx);

    if (nSamples < 2)
    {
        mexWarnMsgIdAndTxt("hnf_mi_mex:nSamples",
                           "fewer than 2 samples; MI will be ill-defined.");
    }

    if (totalElems % nSamples != 0)
    {
        mexErrMsgIdAndTxt("hnf_mi_mex:internal",
                          "Number of elements not divisible by nSamples.");
    }

    nVars = totalElems / nSamples;

    /* Pointers to input data */
    phas = mxGetPr(phas_mx);
    pow  = mxGetPr(pow_mx);

    /* Allocate work arrays */
    sumAmp   = (double *)mxCalloc((size_t)(nbin_i * nVars), sizeof(double));
    countAmp = (mwSize *)mxCalloc((size_t)(nbin_i * nVars), sizeof(mwSize));
    colSum   = (double *)mxCalloc((size_t)nVars, sizeof(double));
    distKL   = (double *)mxCalloc((size_t)nVars, sizeof(double));
    mi_vec   = (double *)mxCalloc((size_t)nVars, sizeof(double));

    if (sumAmp == NULL ||
        countAmp == NULL ||
        colSum == NULL ||
        distKL == NULL ||
        mi_vec == NULL)
    {
        mexErrMsgIdAndTxt("hnf_mi_mex:alloc",
                          "Memory allocation failed.");
    }

    binWidth = (2.0 * PI) / (double)nbin_i;
    log_nbin = log((double)nbin_i);
    NaNval   = mxGetNaN();

    /* --- Main accumulation: binning and amplitude sums --- */

    for (v = 0; v < nVars; ++v)
    {
        mwSize base = v * nSamples;

        for (s = 0; s < nSamples; ++s)
        {
            mwSize idx = base + s;

            double phi   = phas[idx];
            double powval = pow[idx];

            /* Wrap phase to [-pi, pi) via atan2(sin, cos) */
            double sinphi = sin(phi);
            double cosphi = cos(phi);
            double wrapped = atan2(sinphi, cosphi); /* [-pi, pi] */

            /* Shift to [0, 2pi) and bin */
            double shifted = wrapped + PI;
            double binPos  = shifted / binWidth;
            mwSize bIdx    = (mwSize)floor(binPos);

            if ((double)bIdx >= (double)nbin_i)
            {
                bIdx = nbin_i - 1;
            }

            /* linear index for (bin, var) */
            {
                mwSize lin = bIdx * nVars + v;

                sumAmp[lin]   += powval;
                countAmp[lin] += 1;
            }
        }
    }

    /* --- Compute mean amplitude per bin (amplBin) and colSum --- */

    for (v = 0; v < nVars; ++v)
    {
        double csum = 0.0;

        for (b = 0; b < nbin_i; ++b)
        {
            mwSize lin  = b * nVars + v;
            mwSize cnt  = countAmp[lin];
            double ampl;

            if (cnt == 0)
            {
                ampl = 0.0;
            }
            else
            {
                ampl = sumAmp[lin] / (double)cnt;
            }

            sumAmp[lin] = ampl;   /* reuse sumAmp as amplBin */
            csum       += ampl;
        }

        colSum[v] = csum;
    }

    /* --- Compute KL divergence and MI per variable --- */

    for (v = 0; v < nVars; ++v)
    {
        double csum = colSum[v];
        double dKL  = 0.0;

        if (csum == 0.0)
        {
            /* zeroCols: later we set MI = NaN */
            distKL[v] = 0.0;
            continue;
        }

        for (b = 0; b < nbin_i; ++b)
        {
            mwSize lin   = b * nVars + v;
            double ampl  = sumAmp[lin];
            double P     = ampl / csum;

            if (!mxIsFinite(P) || P <= 0.0)
            {
                continue;
            }

            dKL += P * log(P * (double)nbin_i);
        }

        distKL[v] = dKL;
    }

    for (v = 0; v < nVars; ++v)
    {
        double csum = colSum[v];

        if (csum == 0.0)
        {
            mi_vec[v] = NaNval;   /* zeroCols → NaN */
        }
        else
        {
            mi_vec[v] = distKL[v] / log_nbin;
        }
    }

    /* --- Create output array with first dim collapsed to 1 --- */

    dimsOut = (mwSize *)mxCalloc((size_t)nDims, sizeof(mwSize));
    if (dimsOut == NULL)
    {
        mexErrMsgIdAndTxt("hnf_mi_mex:alloc",
                          "Output dimension allocation failed.");
    }

    dimsOut[0] = 1;
    for (i = 1; i < nDims; ++i)
    {
        dimsOut[i] = dims[i];
    }

    mi_mx = mxCreateNumericArray(nDims, dimsOut, mxDOUBLE_CLASS, mxREAL);
    if (mi_mx == NULL)
    {
        mexErrMsgIdAndTxt("hnf_mi_mex:alloc",
                          "Output array allocation failed.");
    }

    mi_out = mxGetPr(mi_mx);

    /* linear layout of mi_vec already matches reshape semantics */
    for (v = 0; v < nVars; ++v)
    {
        mi_out[v] = mi_vec[v];
    }

    plhs[0] = mi_mx;

    /* --- Cleanup --- */

    mxFree(sumAmp);
    mxFree(countAmp);
    mxFree(colSum);
    mxFree(distKL);
    mxFree(mi_vec);
    mxFree(dimsOut);
}