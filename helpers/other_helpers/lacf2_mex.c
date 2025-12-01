#include "mex.h"

/*
 * lacf2_mex(x, mlag)
 *
 * x    : vector (real or complex), length N
 * mlag : maximum lag (scalar, <= N). If omitted or empty, defaults to N.
 *
 * Output:
 *   lacf : mlag x N matrix
 *          lacf(tau, t) = conj(x(t)) * x(t + tau - 1), tau = 1..mtau
 *          (zeros for tau > mtau), exactly as in the original MATLAB code.
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* argument checks */
    if (nrhs < 1 || nrhs > 2)
    {
        mexErrMsgIdAndTxt("lacf2_mex:nrhs",
                          "Usage: lacf = lacf2_mex(x, mlag)");
    }

    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("lacf2_mex:nlhs",
                          "One output argument expected.");
    }

    const mxArray *x_in = prhs[0];

    if (!mxIsDouble(x_in))
    {
        mexErrMsgIdAndTxt("lacf2_mex:type",
                          "Input x must be double.");
    }

    /* flatten x (MATLAB already stores as linear memory) */
    mwSize N = mxGetNumberOfElements(x_in);

    if (N == 0)
    {
        mexErrMsgIdAndTxt("lacf2_mex:empty",
                          "Input x must be non-empty.");
    }

    mwSize mlag;

    if (nrhs >= 2 && !mxIsEmpty(prhs[1]))
    {
        if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
            mxGetNumberOfElements(prhs[1]) != 1)
        {
            mexErrMsgIdAndTxt("lacf2_mex:mlag",
                              "mlag must be a real double scalar.");
        }

        double mlag_d = mxGetScalar(prhs[1]);

        if (mlag_d < 1)
        {
            mexErrMsgIdAndTxt("lacf2_mex:mlag",
                              "mlag must be >= 1.");
        }

        mlag = (mwSize)mlag_d;
    }
    else
    {
        mlag = N;
    }

    if (mlag > N)
    {
        mexErrMsgIdAndTxt("lacf2_mex:mlag",
                          "mlag must be <= length(x).");
    }

    /* detect real vs complex input */
    bool xIsComplex = mxIsComplex(x_in);

    const double *xpr = mxGetPr(x_in);
    const double *xpi = mxGetPi(x_in);

    /* allocate output: mlag x N, real or complex depending on x */
    mxArray *lacf_out;

    if (xIsComplex)
    {
        lacf_out = mxCreateDoubleMatrix(mlag, N, mxCOMPLEX);
    }
    else
    {
        lacf_out = mxCreateDoubleMatrix(mlag, N, mxREAL);
    }

    double *lpr = mxGetPr(lacf_out);
    double *lpi = mxGetPi(lacf_out);

    /* core loop: for t = 1:N
     *              mtau = min(mlag, N - t + 1)
     *              lacf(1:mtau, t) = conj(x(t)) * x(t:t+mtau-1)
     *
     * In zero-based indexing:
     *   t_idx = t - 1, tau_idx = tau - 1
     *   mtau = min(mlag, N - t_idx)
     */

    for (mwSize t_idx = 0; t_idx < N; ++t_idx)
    {
        mwSize mtau = mlag;

        if (N - t_idx < mtau)
        {
            mtau = N - t_idx;
        }

        if (xIsComplex)
        {
            double a = xpr[t_idx];
            double b = xpi[t_idx];

            for (mwSize tau_idx = 0; tau_idx < mtau; ++tau_idx)
            {
                mwSize x2_idx = t_idx + tau_idx;

                double c = xpr[x2_idx];
                double d = xpi[x2_idx];

                double real_part = a * c + b * d;
                double imag_part = a * d - b * c;

                mwSize out_idx = tau_idx + t_idx * mlag;

                lpr[out_idx] = real_part;
                lpi[out_idx] = imag_part;
            }
        }
        else
        {
            double a = xpr[t_idx];

            for (mwSize tau_idx = 0; tau_idx < mtau; ++tau_idx)
            {
                mwSize x2_idx = t_idx + tau_idx;

                double c = xpr[x2_idx];

                mwSize out_idx = tau_idx + t_idx * mlag;

                lpr[out_idx] = a * c;
            }
        }
    }

    plhs[0] = lacf_out;
}