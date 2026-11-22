/* Fixed, zero-peak fits for nf_tfspecparam. Each column is one log10
 * spectrum. With x = -log10(f), OLS gives exponent = cov(x,y)/var(x) and
 * offset = mean(y)-exponent*mean(x). Stream one column at a time in double
 * precision, retain the input precision for surfaces, and summarize the
 * residual only when needed. Fitting, data/fit output and invalid-column
 * semantics mirror nf_fixed_arrays_matlab_local.
 */
#include "mex.h"
#include <math.h>

static int fullRealDouble(const mxArray *a)
{
    return mxIsDouble(a) && !mxIsSparse(a) && !mxIsComplex(a);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize bins, slices, col, f;
    const mxArray *spectra;
    const double *x, *centeredX, *dataDouble = NULL;
    const float *dataSingle = NULL;
    double *apDouble = NULL, *oscDouble = NULL, *summary = NULL;
    float *apSingle = NULL, *oscSingle = NULL;
    mxLogical *badSlice;
    double meanX, sumSquaresX, nanValue = mxGetNaN();
    int returnFit, isSingle, i;

    if (nrhs != 6 || (nlhs != 3 && nlhs != 4))
        mexErrMsgIdAndTxt("NeuroFreq:fixedSpectra:arity","Expected six inputs and three or four outputs.");
    spectra = prhs[0];
    isSingle = mxIsSingle(spectra);
    if ((!isSingle && !mxIsDouble(spectra)) || mxIsComplex(spectra) ||
        mxIsSparse(spectra) || mxGetNumberOfDimensions(spectra) != 2)
        mexErrMsgIdAndTxt("NeuroFreq:fixedSpectra:data","Expected a full real single or double spectra matrix.");
    for (i = 1; i < 6; ++i)
        if (!fullRealDouble(prhs[i]))
            mexErrMsgIdAndTxt("NeuroFreq:fixedSpectra:type","Predictor and options must be full real double arrays.");
    bins = mxGetM(spectra);
    slices = mxGetN(spectra);
    if (bins < 2 || mxGetNumberOfElements(prhs[1]) != (size_t)bins ||
        mxGetNumberOfElements(prhs[2]) != (size_t)bins || mxGetNumberOfElements(prhs[3]) != 1 ||
        mxGetNumberOfElements(prhs[4]) != 1 || mxGetNumberOfElements(prhs[5]) != 1)
        mexErrMsgIdAndTxt("NeuroFreq:fixedSpectra:size","Check predictor length and scalar options.");
    x = mxGetPr(prhs[1]);
    centeredX = mxGetPr(prhs[2]);
    meanX = mxGetScalar(prhs[3]);
    sumSquaresX = mxGetScalar(prhs[4]);
    if (!mxIsFinite(meanX) || !mxIsFinite(sumSquaresX) || sumSquaresX <= 0 ||
        (mxGetScalar(prhs[5]) != 0 && mxGetScalar(prhs[5]) != 1))
        mexErrMsgIdAndTxt("NeuroFreq:fixedSpectra:predictor","Invalid predictor statistics or fit flag.");
    returnFit = mxGetScalar(prhs[5]) != 0;
    plhs[0] = mxCreateNumericMatrix(bins, slices, mxGetClassID(spectra), mxREAL);
    plhs[1] = mxCreateNumericMatrix(bins, slices, mxGetClassID(spectra), mxREAL);
    plhs[2] = mxCreateLogicalMatrix(1, slices);
    badSlice = (mxLogical *)mxGetData(plhs[2]);
    if (nlhs == 4) {
        /* Four contiguous rows per column: offset, exponent, MAE, r squared. */
        plhs[3] = mxCreateDoubleMatrix(4, slices, mxREAL);
        summary = mxGetPr(plhs[3]);
    }
    if (isSingle) {
        dataSingle = (const float *)mxGetData(spectra);
        apSingle = (float *)mxGetData(plhs[0]);
        oscSingle = (float *)mxGetData(plhs[1]);
    } else {
        dataDouble = mxGetPr(spectra);
        apDouble = mxGetPr(plhs[0]);
        oscDouble = mxGetPr(plhs[1]);
    }
    for (col = 0; col < slices; ++col) {
        mwSize start = col*bins;
        double meanY = 0.0, cross = 0.0, sumSquaresY = 0.0, absError = 0.0;
        double exponent, offset;
        for (f = 0; f < bins; ++f) {
            double y = isSingle ? dataSingle[start+f] : dataDouble[start+f];
            if (!isfinite(y)) { badSlice[col] = 1; break; }
            meanY += y;
        }
        if (badSlice[col]) {
            for (f = 0; f < bins; ++f) {
                if (isSingle) apSingle[start+f] = oscSingle[start+f] = (float)nanValue;
                else apDouble[start+f] = oscDouble[start+f] = nanValue;
            }
            if (summary) for (i = 0; i < 4; ++i) summary[4*col+i] = nanValue;
            continue;
        }

        meanY /= (double)bins;
        for (f = 0; f < bins; ++f) {
            double y = isSingle ? dataSingle[start+f] : dataDouble[start+f];
            double centeredY = y - meanY;
            cross += centeredX[f]*centeredY;
            if (summary) sumSquaresY += centeredY*centeredY;
        }
        exponent = cross / sumSquaresX;
        offset = meanY - exponent*meanX;
        for (f = 0; f < bins; ++f) {
            double y = isSingle ? dataSingle[start+f] : dataDouble[start+f];
            double fitted = x[f]*exponent + offset;
            double residual = y - fitted;
            double ap = returnFit ? fitted : y;
            double osc = returnFit ? 0.0 : residual;
            if (isSingle) {
                apSingle[start+f] = (float)ap;
                oscSingle[start+f] = (float)osc;
            } else {
                apDouble[start+f] = ap;
                oscDouble[start+f] = osc;
            }
            if (summary) absError += fabs(residual);
        }
        if (summary) {
            double denominator = sumSquaresY * (exponent*exponent*sumSquaresX);
            double r = denominator > 0 ? (exponent*cross) / sqrt(denominator) : nanValue;
            summary[4*col] = offset;
            summary[4*col+1] = exponent;
            summary[4*col+2] = absError / (double)bins;
            summary[4*col+3] = r*r;
        }
    }
}
