#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    const mxArray *alpha_mx;
    const mxArray *x_mx;
    double *alpha;
    double *x;
    mwSize nDims;
    const mwSize *dims;
    mwSize i;
    mwSize nSamples;
    mwSize nVars;
    mwSize totalElems;
    mwSize v;
    mwSize s;

    double *S;
    double *C;
    double *mx;
    double *ms;
    double *mc;
    double *rxs;
    double *rxc;
    double *rcs;
    double *rho_vec;

    mxArray *rho_mx;
    mxArray *pval_mx;
    double *rho_out;
    double *pval_out;
    mwSize *dimsOut;

    double NaNval;

    if (nrhs != 2)
    {
        mexErrMsgIdAndTxt("hnf_circ_corrcl_mex:nrhs",
                          "Two input arguments required: alpha, x.");
    }
    if (nlhs > 2)
    {
        mexErrMsgIdAndTxt("hnf_circ_corrcl_mex:nlhs",
                          "Too many output arguments.");
    }

    alpha_mx = prhs[0];
    x_mx     = prhs[1];

    if (!mxIsDouble(alpha_mx) || mxIsComplex(alpha_mx))
    {
        mexErrMsgIdAndTxt("hnf_circ_corrcl_mex:alphaType",
                          "alpha must be real double.");
    }
    if (!mxIsDouble(x_mx) || mxIsComplex(x_mx))
    {
        mexErrMsgIdAndTxt("hnf_circ_corrcl_mex:xType",
                          "x must be real double.");
    }

    nDims = mxGetNumberOfDimensions(alpha_mx);
    dims  = mxGetDimensions(alpha_mx);

    if (nDims < 1)
    {
        mexErrMsgIdAndTxt("hnf_circ_corrcl_mex:dimensions",
                          "alpha must have at least one dimension.");
    }

    {
        mwSize nDimsX = mxGetNumberOfDimensions(x_mx);
        const mwSize *dimsX = mxGetDimensions(x_mx);
        if (nDimsX != nDims)
        {
            mexErrMsgIdAndTxt("hnf_circ_corrcl_mex:sizeMismatch",
                              "alpha and x must have the same size.");
        }
        for (i = 0; i < nDims; ++i)
        {
            if (dims[i] != dimsX[i])
            {
                mexErrMsgIdAndTxt("hnf_circ_corrcl_mex:sizeMismatch",
                                  "alpha and x must have the same size.");
            }
        }
    }

    nSamples   = dims[0];
    totalElems = mxGetNumberOfElements(alpha_mx);

    if (nSamples < 2)
    {
        mexErrMsgIdAndTxt("hnf_circ_corrcl_mex:nSamples",
                          "Not enough samples along first dimension.");
    }

    if (totalElems % nSamples != 0)
    {
        mexErrMsgIdAndTxt("hnf_circ_corrcl_mex:internal",
                          "Number of elements not divisible by nSamples.");
    }

    nVars = totalElems / nSamples;

    alpha = mxGetPr(alpha_mx);
    x     = mxGetPr(x_mx);

    S       = (double *)mxCalloc((size_t)totalElems, sizeof(double));
    C       = (double *)mxCalloc((size_t)totalElems, sizeof(double));
    mx      = (double *)mxCalloc((size_t)nVars, sizeof(double));
    ms      = (double *)mxCalloc((size_t)nVars, sizeof(double));
    mc      = (double *)mxCalloc((size_t)nVars, sizeof(double));
    rxs     = (double *)mxCalloc((size_t)nVars, sizeof(double));
    rxc     = (double *)mxCalloc((size_t)nVars, sizeof(double));
    rcs     = (double *)mxCalloc((size_t)nVars, sizeof(double));
    rho_vec = (double *)mxCalloc((size_t)nVars, sizeof(double));

    if (S == NULL ||
        C == NULL ||
        mx == NULL ||
        ms == NULL ||
        mc == NULL ||
        rxs == NULL ||
        rxc == NULL ||
        rcs == NULL ||
        rho_vec == NULL)
    {
        mexErrMsgIdAndTxt("hnf_circ_corrcl_mex:alloc",
                          "Memory allocation failed.");
    }

    for (i = 0; i < totalElems; ++i)
    {
        S[i] = sin(alpha[i]);
        C[i] = cos(alpha[i]);
    }

    for (v = 0; v < nVars; ++v)
    {
        double sumX  = 0.0;
        double sumS  = 0.0;
        double sumC  = 0.0;
        mwSize base  = v * nSamples;

        for (s = 0; s < nSamples; ++s)
        {
            mwSize idx = base + s;
            sumX += x[idx];
            sumS += S[idx];
            sumC += C[idx];
        }

        mx[v] = sumX / (double)nSamples;
        ms[v] = sumS / (double)nSamples;
        mc[v] = sumC / (double)nSamples;
    }

    NaNval = mxGetNaN();

    for (v = 0; v < nVars; ++v)
    {
        double sumXcSc = 0.0;
        double sumXc2  = 0.0;
        double sumSc2  = 0.0;

        double sumXcCc = 0.0;
        double sumCc2  = 0.0;
        double sumCcSc = 0.0;

        double mx_v    = mx[v];
        double ms_v    = ms[v];
        double mc_v    = mc[v];

        mwSize base    = v * nSamples;

        for (s = 0; s < nSamples; ++s)
        {
            mwSize idx = base + s;

            double Xval = x[idx];
            double Sval = S[idx];
            double Cval = C[idx];

            double Xc = Xval - mx_v;
            double Sc = Sval - ms_v;
            double Cc = Cval - mc_v;

            sumXcSc += Xc * Sc;
            sumXc2  += Xc * Xc;
            sumSc2  += Sc * Sc;

            sumXcCc += Xc * Cc;
            sumCc2  += Cc * Cc;
            sumCcSc += Cc * Sc;
        }

        {
            double den_xs = sqrt(sumXc2 * sumSc2);
            double den_xc = sqrt(sumXc2 * sumCc2);
            double den_cs = sqrt(sumCc2 * sumSc2);

            if (den_xs == 0.0)
            {
                rxs[v] = NaNval;
            }
            else
            {
                rxs[v] = sumXcSc / den_xs;
            }

            if (den_xc == 0.0)
            {
                rxc[v] = NaNval;
            }
            else
            {
                rxc[v] = sumXcCc / den_xc;
            }

            if (den_cs == 0.0)
            {
                rcs[v] = NaNval;
            }
            else
            {
                rcs[v] = sumCcSc / den_cs;
            }
        }
    }

    for (v = 0; v < nVars; ++v)
    {
        double rxs_v = rxs[v];
        double rxc_v = rxc[v];
        double rcs_v = rcs[v];

        double den_rho;
        double num_rho;
        double frac;

        den_rho = 1.0 - rcs_v * rcs_v;

        if (den_rho <= 0.0)
        {
            rho_vec[v] = NaNval;
        }
        else
        {
            num_rho = (rxc_v * rxc_v) +
                      (rxs_v * rxs_v) -
                      2.0 * rxc_v * rxs_v * rcs_v;

            frac = num_rho / den_rho;

            if (frac < 0.0)
            {
                frac = 0.0;
            }

            rho_vec[v] = sqrt(frac);
        }
    }

    dimsOut = (mwSize *)mxCalloc((size_t)nDims, sizeof(mwSize));
    if (dimsOut == NULL)
    {
        mexErrMsgIdAndTxt("hnf_circ_corrcl_mex:alloc",
                          "Memory allocation failed.");
    }

    dimsOut[0] = 1;
    for (i = 1; i < nDims; ++i)
    {
        dimsOut[i] = dims[i];
    }

    rho_mx  = mxCreateNumericArray(nDims, dimsOut, mxDOUBLE_CLASS, mxREAL);
    pval_mx = mxCreateNumericArray(nDims, dimsOut, mxDOUBLE_CLASS, mxREAL);

    if (rho_mx == NULL || pval_mx == NULL)
    {
        mexErrMsgIdAndTxt("hnf_circ_corrcl_mex:alloc",
                          "Output allocation failed.");
    }

    rho_out  = mxGetPr(rho_mx);
    pval_out = mxGetPr(pval_mx);

    for (v = 0; v < nVars; ++v)
    {
        double rho_v = rho_vec[v];
        double p_v;

        rho_out[v] = rho_v;

        if (rho_v != rho_v)
        {
            p_v = NaNval;
        }
        else
        {
            double z = ((double)nSamples) * rho_v * rho_v;
            p_v = exp(-0.5 * z);
        }

        pval_out[v] = p_v;
    }

    plhs[0] = rho_mx;
    if (nlhs > 1)
    {
        plhs[1] = pval_mx;
    }
    else
    {
        mxDestroyArray(pval_mx);
    }

    mxFree(S);
    mxFree(C);
    mxFree(mx);
    mxFree(ms);
    mxFree(mc);
    mxFree(rxs);
    mxFree(rxc);
    mxFree(rcs);
    mxFree(rho_vec);
    mxFree(dimsOut);
}