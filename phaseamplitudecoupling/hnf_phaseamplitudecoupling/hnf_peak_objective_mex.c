/* Scalar peak-fit error and optional exact gradient for hnf_fit_peak_guess.
 * A column-major N-by-3 parameter matrix (or its packed vector) contains the
 * centers, heights and widths. Keeping just one bin's N basis values avoids
 * the peak-by-frequency arrays otherwise allocated on every solver evaluation.
 *
 * Inputs: parameters, frequencies, observations, starting guess, penalty
 * weight, shape (0 = Gaussian, 1 = Cauchy). A zero weight ignores the guess.
 * The Matlab fallback implements the same formulas in hnf_error_model and
 * hnf_error_model_constr. The optimizer and its stopping conditions stay there.
 */
#include "mex.h"
#include <math.h>
#include <float.h>

static int fullRealDouble(const mxArray *a)
{
    return mxIsDouble(a) && !mxIsComplex(a) && !mxIsSparse(a);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t n, bins, k, j;
    const double *params, *freqs, *observed, *guess;
    double weight, shape, error = 0.0;
    double *gradient = NULL, *basis = NULL, *values = NULL;
    int i;

    if (nrhs != 6 || nlhs < 1 || nlhs > 2)
        mexErrMsgIdAndTxt("NeuroFreq:peakObjective:arity",
            "Expected six inputs and one or two outputs.");
    for (i = 0; i < nrhs; ++i)
        if (!fullRealDouble(prhs[i]))
            mexErrMsgIdAndTxt("NeuroFreq:peakObjective:type",
                "All inputs must be full real double arrays.");

    bins = mxGetNumberOfElements(prhs[1]);
    n = mxGetNumberOfElements(prhs[0]) / 3;
    if (mxGetNumberOfElements(prhs[0]) != 3*n ||
        mxGetNumberOfElements(prhs[2]) != bins ||
        mxGetNumberOfElements(prhs[4]) != 1 ||
        mxGetNumberOfElements(prhs[5]) != 1)
        mexErrMsgIdAndTxt("NeuroFreq:peakObjective:size",
            "Expected parameter triples, equally long spectra and scalar options.");
    weight = mxGetScalar(prhs[4]);
    shape = mxGetScalar(prhs[5]);
    if (!(weight >= 0.0) || !mxIsFinite(weight) || (shape != 0.0 && shape != 1.0) ||
        (weight > 0.0 && mxGetNumberOfElements(prhs[3]) != 3*n))
        mexErrMsgIdAndTxt("NeuroFreq:peakObjective:options",
            "Check the penalty weight, starting guess dimensions and peak shape.");

    params = mxGetPr(prhs[0]);
    freqs = mxGetPr(prhs[1]);
    observed = mxGetPr(prhs[2]);
    guess = mxGetPr(prhs[3]);
    if (nlhs == 2) {
        plhs[1] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),
            mxGetDimensions(prhs[0]), mxDOUBLE_CLASS, mxREAL);
        gradient = mxGetPr(plhs[1]);
        basis = (double *)mxCalloc(n, sizeof(double));
        values = (double *)mxCalloc(n, sizeof(double));
    }

    for (j = 0; j < bins; ++j) {
        double model = 0.0, residual;
        for (k = 0; k < n; ++k) {
            double u = (freqs[j] - params[k]) / params[2*n + k];
            double b = shape == 1.0 ? 1.0 / (1.0 + u*u) : exp(-0.5*u*u);
            double height = params[n+k], peak = height*b;
            if (shape == 1.0 && fabs(u) > sqrt(DBL_MAX)/2.0) {
                double inv = params[2*n+k]/(freqs[j]-params[k]);
                b = (inv*inv)/(1.0+inv*inv);
                peak = ((height*inv)*inv)/(1.0+inv*inv);
            } else if (shape==0.0 && fabs(height)>1.0 && -0.5*u*u < log(DBL_MIN)) {
                peak = copysign(exp(log(fabs(height))-0.5*u*u),height);
            }
            model += peak;
            if (basis) { basis[k] = b; values[k] = peak; }
        }
        residual = observed[j] - model;
        error += residual*residual;
        if (gradient) {
            for (k = 0; k < n; ++k) {
                double width = params[2*n + k];
                double u = (freqs[j] - params[k])/width;
                double peak = values[k], center, scale;
                if (shape == 1.0) {
                    double q, w;
                    if (fabs(u)>1.0) {
                        double inv = 1.0/u;
                        w = 1.0/(1.0+inv*inv);
                        q = inv*w;
                    } else {
                        q = u*basis[k];
                        w = u*q;
                    }
                    center = (2.0*(peak*q))/width;
                    scale = (2.0*(peak*w))/width;
                } else {
                    if (!isfinite(u) && basis[k] == 0) u = 0;
                    center = (peak*u)/width;
                    scale = ((peak*u)*u)/width;
                }
                gradient[k]       -= 2.0*residual*center;
                gradient[n + k]   -= 2.0*residual*basis[k];
                gradient[2*n + k] -= 2.0*residual*scale;
            }
        }
    }

    if (weight > 0.0) {
        /* Center and height have separate weak/strong penalties. Width is free. */
        for (j = 0; j < 2; ++j) {
            double sumSquares = 0.0;
            for (k = j*n; k < (j+1)*n; ++k) {
                double d = params[k] - guess[k];
                sumSquares += d*d;
                if (gradient) gradient[k] += 2.0*weight*d;
            }
            error += weight*sumSquares;
        }
    }
    if (basis) { mxFree(basis); mxFree(values); }
    plhs[0] = mxCreateDoubleScalar(error);
}
