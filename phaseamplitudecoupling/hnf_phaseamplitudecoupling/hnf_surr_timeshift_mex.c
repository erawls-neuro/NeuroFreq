#include "mex.h"
#include <math.h>

/*
 * Xs = hnf_surr_timeshift_mex(Xmat, shifts)
 *
 * INPUT:
 *   Xmat   : double matrix [nSamples x nVars]
 *   shifts : double vector [1 x nVars] or [nVars x 1]
 *
 * OUTPUT:
 *   Xs     : double matrix [nSamples x nVars],
 *            where each column v is circularly shifted by shifts(v)
 *            samples along the first dimension.
 *
 * Notes:
 *   - Shifts are assumed to be integers in [0, nSamples-1].
 *   - Data is assumed column-major (standard MATLAB).
 */

void mexFunction( int nlhs,
                  mxArray *plhs[],
                  int nrhs,
                  const mxArray *prhs[] )
{
    /* Check number of inputs and outputs */
    if( nrhs != 2 )
    {
        mexErrMsgIdAndTxt( "hnf_surr_timeshift_mex:nrhs",
                           "Two input arguments required: Xmat, shifts." );
    }

    if( nlhs > 1 )
    {
        mexErrMsgIdAndTxt( "hnf_surr_timeshift_mex:nlhs",
                           "One output argument allowed." );
    }

    /* Get Xmat */
    const mxArray *Xmat_mx = prhs[0];
    if( !mxIsDouble(Xmat_mx) || mxIsComplex(Xmat_mx) )
    {
        mexErrMsgIdAndTxt( "hnf_surr_timeshift_mex:Xmat",
                           "Xmat must be a real double matrix." );
    }

    mwSize nSamples = mxGetM(Xmat_mx);
    mwSize nVars    = mxGetN(Xmat_mx);

    if( nSamples == 0 || nVars == 0 )
    {
        /* Degenerate case: just return a copy */
        plhs[0] = mxDuplicateArray(Xmat_mx);
        return;
    }

    /* Get shifts */
    const mxArray *shifts_mx = prhs[1];
    if( !mxIsDouble(shifts_mx) || mxIsComplex(shifts_mx) )
    {
        mexErrMsgIdAndTxt( "hnf_surr_timeshift_mex:shifts",
                           "shifts must be a real double vector." );
    }

    mwSize sM = mxGetM(shifts_mx);
    mwSize sN = mxGetN(shifts_mx);
    mwSize nShifts = sM * sN;

    if( nShifts != nVars )
    {
        mexErrMsgIdAndTxt( "hnf_surr_timeshift_mex:shiftsSize",
                           "Number of shifts must match number of columns in Xmat." );
    }

    const double *Xmat   = mxGetPr(Xmat_mx);
    const double *shifts = mxGetPr(shifts_mx);

    /* Create output array */
    plhs[0] = mxCreateDoubleMatrix(nSamples, nVars, mxREAL);
    double *Xs = mxGetPr(plhs[0]);

    /* Column-major circular shift */
    for( mwSize v = 0; v < nVars; v++ )
    {
        double shift_val = shifts[v];

        /* Round to nearest integer, since MATLAB gave us double */
        long k_long = (long) llround(shift_val);

        if( k_long < 0 || k_long >= (long)nSamples )
        {
            mexErrMsgIdAndTxt( "hnf_surr_timeshift_mex:shiftRange",
                               "Each shift must be in [0, nSamples-1]." );
        }

        mwSize k = (mwSize) k_long;

        const double *x_col  = Xmat + v * nSamples;
        double       *xs_col = Xs   + v * nSamples;

        if( k == 0 )
        {
            /* No shift needed, copy column directly */
            for( mwSize r = 0; r < nSamples; r++ )
            {
                xs_col[r] = x_col[r];
            }
        }
        else
        {
            /* Circular shift down by k samples:
             * xs(r) = x( (r - k) mod nSamples )
             *
             * We implement as:
             *   src = r - k
             *   if src < 0: src += nSamples
             */
            for( mwSize r = 0; r < nSamples; r++ )
            {
                long src = (long)r - (long)k;
                if( src < 0 )
                {
                    src += (long)nSamples;
                }
                xs_col[r] = x_col[src];
            }
        }
    }
}