/* Finite EEG window, channels x samples. Return the sample-wise common-average
 * reference and the three QC statistics needed for artifact classification.
 * The caller keeps NaN-bearing windows on the omission-aware MATLAB path. */
#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <float.h>

static int ordered(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x>y)-(x<y);
}

static void add(double x, double *sum, double *correction)
{
    double next = *sum + x;
    *correction += (fabs(*sum) >= fabs(x)) ? (*sum-next)+x : (x-next)+*sum;
    *sum = next;
}

void mexFunction(int nlhs, mxArray *out[], int nrhs, const mxArray *arg[])
{
    if (nrhs!=1 || nlhs!=4 || !mxIsDouble(arg[0]) || mxIsSparse(arg[0]) ||
        mxIsComplex(arg[0]) || mxGetNumberOfDimensions(arg[0])!=2 || mxIsEmpty(arg[0]))
        mexErrMsgIdAndTxt("NeuroFreq:QCWindow", "Expected a nonempty full real double matrix and four outputs.");
    const mwSize channels = mxGetM(arg[0]), samples = mxGetN(arg[0]);
    const double *x = mxGetPr(arg[0]);
    out[0] = mxCreateDoubleMatrix(channels,samples,mxREAL);
    for (int j=1; j<4; ++j) out[j] = mxCreateDoubleMatrix(channels,1,mxREAL);
    double *ref = mxGetPr(out[0]), *sd = mxGetPr(out[1]);
    double *peak = mxGetPr(out[2]), *rms = mxGetPr(out[3]);
    double *row = mxCalloc(samples,sizeof(double));
    for (mwSize t=0; t<samples; ++t) {
        const double *column = x+t*channels;
        double scale=0, sum=0, correction=0;
        for (mwSize c=0; c<channels; ++c) {
            if (!mxIsFinite(column[c]) || fabs(column[c])>DBL_MAX/4)
                mexErrMsgIdAndTxt("NeuroFreq:QCWindowRange","Finite native windows require |x| <= realmax/4.");
            if (fabs(column[c])>scale) scale=fabs(column[c]);
        }
        if (scale==0) continue;
        for (mwSize c=0; c<channels; ++c) add(column[c]/scale,&sum,&correction);
        double mean=((sum+correction)/channels)*scale;
        for (mwSize c=0; c<channels; ++c) ref[c+t*channels]=column[c]-mean;
    }
    for (mwSize c=0; c<channels; ++c) {
        double sum=0, correction=0, scale=0, origin=x[c];
        for (mwSize t=0; t<samples; ++t) {
            double deviation=fabs(x[c+t*channels]-origin);
            if (deviation>scale) scale=deviation;
            row[t]=ref[c+t*channels];
        }
        if (scale!=0 && samples>1) {
            for (mwSize t=0; t<samples; ++t)
                add((x[c+t*channels]-origin)/scale,&sum,&correction);
            double mean=(sum+correction)/samples;
            sum=0; correction=0;
            for (mwSize t=0; t<samples; ++t) {
                double deviation=(x[c+t*channels]-origin)/scale-mean;
                add(deviation*deviation,&sum,&correction);
            }
            sd[c]=sqrt((sum+correction)/(samples-1))*scale;
        }
        qsort(row,samples,sizeof(double),ordered);
        double a=row[(samples-1)/2], b=row[samples/2];
        double median=a+0.5*(b-a);
        peak[c]=fmax(median-row[0],row[samples-1]-median);
        if (peak[c]==0) continue;
        sum=0; correction=0;
        for (mwSize t=0; t<samples; ++t) {
            double v=(row[t]-median)/peak[c];
            add(v*v,&sum,&correction);
        }
        rms[c]=sqrt((sum+correction)/samples)*peak[c];
    }
    mxFree(row);
}
