/* mex_dtspg.c
 * Source code for the mex_dtspg.mex file. This MEX file is intended to be
 * invoked the following way:
 *
 *   x = mex_dtspg(t, b)
 */

#include "mex.h"
#include <csp_dtspg.h>

/* Entry point */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int     n = mxGetN(prhs[0]);
  int  bSize= mxGetM(prhs[0]);
  double *t = mxGetPr(prhs[0]);
  double *x;
  double *b = mxGetPr(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
  x = mxGetPr(plhs[0]); 

  csp_dtspg_sv( n, bSize, t, b, x, mxGetN(prhs[1]));
} 
